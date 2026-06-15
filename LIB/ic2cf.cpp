#include "gaboom.h"
#include "fileio.h"
#include "LigandRingFlex/LigandRingFlex.h"   // Phase 2: ring pucker apply

#ifdef _OPENMP
#  include <omp.h>
#endif

// A3 perf: file-scope struct definitions so thread_local scratch can hold them.
struct IC2CFSavedAtom      { int idx; atom value; };
struct IC2CFSavedResRot    { int idx; int  rot;   };

// Per-thread reusable buffers — eliminates per-eval heap allocs under OMP.
struct IC2CFScratch {
	std::vector<IC2CFSavedAtom>           saved_atoms;      // reserve 64
	std::vector<IC2CFSavedResRot>         saved_res_rots;   // reserve 16
	std::vector<char>                     saved_res_seen;   // sized to res_cnt+1
	std::vector<std::pair<int,int>>       intraclashes;
};
// One IC2CFScratch per OMP thread (no synchronisation needed — each thread
// has exclusive access to its own scratch throughout ic2cf execution).
static thread_local IC2CFScratch tl_ic2cf_scratch;

/******************************************************************************
 * SUBROUTINE ic2cf gets a vector with internal coordinates rebuilds the 
 * cartesian coordinates and calculates the complementarity function. Its 
 * input vector has the list of ic's that are to be optimized and a global
 * vector contains the information of what kind of variable each item in icv
 * is and to which residue it belongs.
 *****************************************************************************/

//THE PROCEDURE SHOULD RECEIVE A 2ND SET OF GENES THAT ENCODES FOR THE ROTAMER DISTRIBUTION IN THE BPK

cfstr ic2cf(FA_Global* FA,VC_Global* VC,atom* atoms,resid* residue,
			gridpoint* cleftgrid,int npar, double* icv)
{
  
	// static int nbranch = 0;
	
	int i,j,k;
	int cat;    /* atom number constrained to the one considered */

	cfstr cf;
	
	int rclash=0;

	psFlexDEE_Node psFlexDEENode;
	int dee_val;
	int rotflag;
  
	unsigned int grd_idx;
	unsigned int rot_idx;
	//float threshold=10.0;
	//int nflxchk=-1;
	int normalmode=-1;
	int deelig_list[100];
	
	//int rigid_clash=0;
	//float rand=0.0;
	//float min_dis=10.0;
  
	// copy values from icv into respective srtructure atom ic fields 
	// andcompute the ic of a constrained atom prior to reconstruction
	
	//for(i=0;i<npar;i++){printf("[%8.3f]",icv[i]);}printf("\n");
	//PAUSE;

	//printf("NEW INDIVIDUAL=");
	for(i=0;i<npar;i++){
		//printf("[%8.3f]",icv[i]);

		if(FA->map_par[i].typ==-1) { //by index
			
			grd_idx = (uint)icv[i];
			atoms[FA->map_par[i].atm].dis = cleftgrid[grd_idx].dis;
			atoms[FA->map_par[i].atm].ang = cleftgrid[grd_idx].ang;
			atoms[FA->map_par[i].atm].dih = cleftgrid[grd_idx].dih;
			
		}else if(FA->map_par[i].typ==0)  {
			atoms[FA->map_par[i].atm].dis = (float)icv[i];
			
		}else if(FA->map_par[i].typ==1)  {
			atoms[FA->map_par[i].atm].ang = (float)icv[i];
			
		}else if(FA->map_par[i].typ==2)  {
			atoms[FA->map_par[i].atm].dih = (float)icv[i];
			
			j=FA->map_par[i].atm;
			cat=atoms[j].rec[3];
			if(cat != 0){
				while(cat != FA->map_par[i].atm){
					atoms[cat].dih=atoms[j].dih + atoms[cat].shift; 
					j=cat;
					cat=atoms[j].rec[3];
				}
			}
			
		}else if(FA->map_par[i].typ==3) { //by index
			grd_idx = (uint)icv[i];
			//printf("icv(index): %d\n", grd_idx);
			//PAUSE;
      
			// serves as flag , but also as grid index
			normalmode=grd_idx;
      
		}else if(FA->map_par[i].typ==4)  {
			rot_idx = (uint)(icv[i]+0.5f);
      
			residue[atoms[FA->map_par[i].atm].ofres].rot=(int)rot_idx;
      
			/*
			  printf("residue[%d].rot[%d] - fatm=%d - latm=%d\n",
			  residue[atoms[FA->map_par[i].atm].ofres].number,
			  residue[atoms[FA->map_par[i].atm].ofres].rot,
			  residue[atoms[FA->map_par[i].atm].ofres].fatm[rot_idx],
			  residue[atoms[FA->map_par[i].atm].ofres].latm[rot_idx]);
			*/
      
		}
    
	}
	//printf("HERE\n");
	//PAUSE;
  
	// do not alter default (ini) protein conf.
	if(normalmode > -1){
		alter_mode(atoms,residue,FA->normal_grid[normalmode],FA->res_cnt,FA->normal_modes);
	}

	// Save FA->ori and the mutable atom/residue state before any modification.
	// ic2cf() is called repeatedly on shared atom/residue buffers in serial
	// paths, so every evaluation must leave the caller's baseline unchanged.
	float ori_save[3] = {FA->ori[0], FA->ori[1], FA->ori[2]};
	// A3: use thread-local pre-allocated scratch (no heap alloc per eval)
	IC2CFScratch& scr = tl_ic2cf_scratch;
	scr.saved_atoms.clear();
	scr.saved_res_rots.clear();
	if(scr.saved_atoms.capacity() < 64)    scr.saved_atoms.reserve(64);
	if(scr.saved_res_rots.capacity() < 16) scr.saved_res_rots.reserve(16);
	// ensure saved_res_seen is large enough and zeroed for this call
	const int res_cnt_p1 = FA->res_cnt + 1;
	if((int)scr.saved_res_seen.size() < res_cnt_p1)
		scr.saved_res_seen.assign(res_cnt_p1, 0);
	else
		std::fill(scr.saved_res_seen.begin(), scr.saved_res_seen.begin() + res_cnt_p1, 0);
	std::vector<IC2CFSavedAtom>&   saved_atoms    = scr.saved_atoms;
	std::vector<IC2CFSavedResRot>& saved_res_rots = scr.saved_res_rots;
	std::vector<char>&             saved_res_seen  = scr.saved_res_seen;
	for (int r = 0; r < FA->nors; ++r)
		for (int m = 0; m < FA->nmov[r]; ++m) {
			int ai = FA->mov[r][m];
			saved_atoms.push_back(IC2CFSavedAtom{ai, atoms[ai]});
		}
	for (i = 0; i < npar; ++i) {
		if (FA->map_par[i].typ == 4) {
			int ri = atoms[FA->map_par[i].atm].ofres;
			if (ri >= 0 && ri <= FA->res_cnt && !saved_res_seen[ri]) {
				saved_res_rots.push_back(IC2CFSavedResRot{ri, residue[ri].rot});
				saved_res_seen[ri] = 1;
			}
		}
	}

	// ── Ring pucker apply (LigandRingFlex Phase 2) ───────────────────────────
	// Snap each furanose ring's internal dihedrals (.dih) to the Cremer-Pople
	// pucker phase carried by the current chromosome (loaded into
	// FA->ring_cur_phases before this call). Must run AFTER the icv→.dih copy
	// loop above and BEFORE buildcc(), so the reconstructed Cartesian coords
	// reflect the pucker. Gated OFF by default; ring bonds are excluded from
	// map_par, so this never perturbs the standard torsional genes.
	if (FA->ring_flex_active && FA->ring_flex_template &&
	    FA->ring_flex_template->has_rings() && FA->ring_n_sugars > 0) {
		const ligand_ring_flex::RingFlexGenes& tmpl = *FA->ring_flex_template;
		std::vector<float> phases(FA->ring_cur_phases,
		                          FA->ring_cur_phases + FA->ring_n_sugars);
		sugar_pucker::apply_sugar_puckers(
			atoms, tmpl.sugar_ring_indices, phases, tmpl.sugar_types);
	}

	/* rebuild cartesian coordinates of optimized residues*/
	for(i=0;i<FA->nors;i++){ //number of optimized residues
		buildcc(FA,atoms,FA->nmov[i],FA->mov[i]);
	}

	// Out-of-bounds penalty: if any moved atom lands >200Å beyond the protein
	// bounding box, the ligand has escaped the grid.  Restore FA->ori and the
	// moved-atom coordinates so serial callers do not inherit corrupted state,
	// then return maximum penalty — the chromosome stays in the population but
	// ranks last, and the GA evolves away from it naturally.
	{
		const float margin = 200.0f;
		bool oob = false;
		for (int r = 0; r < FA->nors && !oob; ++r) {
			for (int m = 0; m < FA->nmov[r] && !oob; ++m) {
				int ai = FA->mov[r][m];
				for (int j = 0; j < 3; ++j) {
					if (atoms[ai].coor[j] < FA->globalmin[j] - margin ||
					    atoms[ai].coor[j] > FA->globalmax[j] + margin) {
						oob = true;
						break;
					}
				}
			}
		}
		if (oob) {
			FA->ori[0] = ori_save[0];
			FA->ori[1] = ori_save[1];
			FA->ori[2] = ori_save[2];
			for (const auto& sa : saved_atoms) {
				atoms[sa.idx] = sa.value;
			}
			for (const auto& sr : saved_res_rots) {
				residue[sr.idx].rot = sr.rot;
			}
			cfstr cf_oob{};
			cf_oob.com = 99999.0;
			return cf_oob;
		}
	}

	// A3: reuse thread-local intraclashes buffer
	scr.intraclashes.clear();
	std::vector<std::pair<int,int>>& intraclashes = scr.intraclashes;
	bool error;
	double penalty = vcfunction(FA,VC,atoms,residue,intraclashes,&error);
	if(error){
		// Fix: named-field init — aggregate order was wrong (metal_coord=1, rclash=0).
		// wal=penalty is the only nonzero energy field; rclash=1 marks it as a clash.
		cfstr cf_clash{};
		cf_clash.wal    = penalty;
		cf_clash.rclash = 1;
		return cf_clash;
	}
	
	cf.com = 0.0;
	cf.wal = 0.0;
	cf.sas = 0.0;
	cf.con = 0.0;
	cf.elec = 0.0;
	cf.hbond = 0.0;
	cf.gist_desolv = 0.0;
	cf.metal_coord = 0.0;
	cf.rclash = 0;
    
	for(i=0;i<FA->num_optres;i++){
    
		resid* res = &residue[FA->optres[i].rnum];
		
		// flexible side-chain optimization
		if ( !FA->optres[i].type ) {
  
			if ( FA->optres[i].cf.rclash == 1 ) { 

				/*
				  printf("%s %c %d is clashing\n",
				  residue[FA->optres[i].rnum].name,
				  residue[FA->optres[i].rnum].chn,
				  residue[FA->optres[i].rnum].number);
				*/
	
				rclash = 1; 
				
			}
      
		}else{
			
			//int fatm = res->fatm[0];
			if(FA->deelig_flex){
				std::vector< std::pair<int,int> >::iterator it;
				for(it=intraclashes.begin(); it!=intraclashes.end(); ++it)
				{
					for(k=1; k<=res->fdih; k++){
						deelig_list[k] = -1000;
					}
					
					// flex bonds list
					int fbindex = 0;
					int* fblist = res->shortflex[it->first][it->second];
					
					//printf("between[%d][%d]\n", atoms[it->first+fatm].number, atoms[it->second+fatm].number);
					//cout << fblist[fbindex] << endl;
					/*
					printf("fblist = [");
					while(fblist[fbindex] != -1){
						printf("%d,", fblist[fbindex]);
						fbindex++;
					}
					printf("]\n");
					*/
					
					fbindex = 0;
					while(fblist[fbindex] != -1){
						if(atoms[res->bond[fblist[fbindex]]].par != NULL){
							deelig_list[fblist[fbindex]] =
								(int)(atoms[res->bond[fblist[fbindex]]].dih + 0.5);
						}
						fbindex++;
					}
										
					struct deelig_node_struct* node = FA->deelig_root_node;
				[[maybe_unused]] bool add = false;
					
						for(k=1; k<=res->fdih; k++){
						std::map<int, struct deelig_node_struct*>::iterator it;
						it = node->childs.find(deelig_list[k]);
						
						if(it == node->childs.end()){
							struct deelig_node_struct* deelig_child_node = new struct deelig_node_struct;
							
							//if(k==1) cout << "new node added " << deelig_list[k] << endl;
							node->childs[deelig_list[k]] = deelig_child_node;
							
							deelig_child_node->parent = node;
							node = deelig_child_node;
							add = true;
						}else{
							node = it->second;
						}
					}
					
					/*
					if(add) { 
						printf("deelig list = [");
						for(k=1; k<=res->fdih; k++){
							printf("%d,", deelig_list[k]);
						}
						printf("]\n");
						
						nbranch++; cout << "total branches " << nbranch << endl;
					}
					*/
				}
			}
		}
		
		/*
		  printf("optres[%2d].cf  .wal = %.3f\n               .com = %.3f\n               .sas = %.3f\n               .con = %.3f\n",
		  i,FA->optres[i].cf.wal,FA->optres[i].cf.com,FA->optres[i].cf.sas,FA->optres[i].cf.con);
		*/
        
		//sum += (FA->optres[i].cf.com - FA->optres[i].cf.wal + FA->optres[i].cf.sas - FA->optres[i].cf.con);
    
		cf.com += FA->optres[i].cf.com;
		cf.wal += FA->optres[i].cf.wal;
		cf.sas += FA->optres[i].cf.sas;
		cf.con += FA->optres[i].cf.con;
		cf.metal_coord += FA->optres[i].cf.metal_coord;
		cf.hbond += FA->optres[i].cf.hbond;

	}

  
	// add rotamer list to dee list
	// When running inside an OpenMP parallel region each thread has its own
	// copy of FA (thread-local), so DEE updates are serialised via a critical
	// section and written back to FA directly (the thread-local FA shares the
	// psFlexDEENode pointer with the master; the critical section prevents
	// concurrent linked-list corruption).
#ifdef _OPENMP
	if (FA->useflexdee > 0 && rclash && !omp_in_parallel()) {
#else
	if (FA->useflexdee > 0 && rclash) {
#endif
    
		NEW( psFlexDEENode, sFlexDEE_Node );

		psFlexDEENode->rotlist = (int*)malloc(FA->nflxsc_real*sizeof(int));
    
    
		// fill rotamer list
		k=0;
		rotflag=0;
    
		for(j=0;j<FA->nflxsc;j++){
      
			if(residue[FA->flex_res[j].inum].trot > 0   &&
			   FA->flex_res[j].cflag != 0){
	
				psFlexDEENode->rotlist[k++] = residue[FA->flex_res[j].inum].rot;
	
				if(residue[FA->flex_res[j].inum].rot != 0) { rotflag=1; }
			}
      
		}
    
    
		// do not add initial conformation to DEE list
		if ( rotflag ) {
      
			/*
			  printf("\n-----------------\nCreating new node...\n");
			  printf("DEE to add = ");for(k=0;k<FA->nflxsc_real;k++){printf("%3d",psFlexDEENode->rotlist[k]);}printf("\n");
	
			  dee_print(FA->psFlexDEENode,FA->nflxsc_real);
	
			  //getchar();
			  */
      
    
			if( FA->psFlexDEENode ) {
	
				//FA->psFlexDEENode = FA->psFlexDEENode->last;
	
				while ( FA->psFlexDEENode->next != NULL ) {
					FA->psFlexDEENode = FA->psFlexDEENode->next;
				}
	
				dee_val = dee_pivot(psFlexDEENode,&FA->psFlexDEENode,1,FA->FlexDEE_Nodes,(int)((FA->FlexDEE_Nodes+1)/2),FA->FlexDEE_Nodes,FA->nflxsc_real);
	
				if ( dee_val == 1 ) {
	  
					if ( FA->psFlexDEENode->next == NULL ) {
	    
						psFlexDEENode->next = NULL;
						psFlexDEENode->prev = FA->psFlexDEENode;
						FA->psFlexDEENode->next = psFlexDEENode;
	    
						psFlexDEENode->first = FA->psFlexDEENode;
	    
						dee_last(FA->psFlexDEENode,psFlexDEENode);
	    
					} else {
	    
						psFlexDEENode->first = FA->psFlexDEENode->first;
						psFlexDEENode->last = FA->psFlexDEENode->last;
	    
						psFlexDEENode->next = FA->psFlexDEENode->next;
						psFlexDEENode->prev = FA->psFlexDEENode;
						FA->psFlexDEENode->next = psFlexDEENode; 
						psFlexDEENode->next->prev = psFlexDEENode;
	    
					}
	  
					FA->FlexDEE_Nodes++;
	  
				} else if ( dee_val == -1 ) {
	  
					if ( FA->psFlexDEENode->prev == NULL ) {
	    
						psFlexDEENode->prev = NULL;
						psFlexDEENode->next = FA->psFlexDEENode;
						FA->psFlexDEENode->prev = psFlexDEENode;
	    
						psFlexDEENode->last = FA->psFlexDEENode;
	    
						dee_first(FA->psFlexDEENode,psFlexDEENode);
	    
					} else {
	    
						psFlexDEENode->first = FA->psFlexDEENode->first;
						psFlexDEENode->last = FA->psFlexDEENode->last;
	    
						psFlexDEENode->prev = FA->psFlexDEENode->prev;
						psFlexDEENode->next = FA->psFlexDEENode;
						FA->psFlexDEENode->prev = psFlexDEENode; 
						psFlexDEENode->prev->next = psFlexDEENode;
	    
					}
	  
					FA->FlexDEE_Nodes++;
	  
				} else {
	  
					FREE(psFlexDEENode);
	  
				}
	
			} else {
	
				FA->psFlexDEENode = psFlexDEENode;
	
				FA->psFlexDEENode->next = NULL;
				FA->psFlexDEENode->prev = NULL;
	
				FA->psFlexDEENode->first = FA->psFlexDEENode;
				FA->psFlexDEENode->last = FA->psFlexDEENode;
	
				FA->FlexDEE_Nodes++;
	
			}
      
		}
    
	}

	// Restore FA->ori only — NOT atoms[] — on the normal scoring exit.
	//
	// FA->ori is global scratch that drifts cumulatively if not reset; always
	// restore it so the next ic2cf call uses the correct receptor-centre frame.
	//
	// atoms[] and residue[].rot are intentionally NOT restored here.
	// The GA parallel path works on thread-private copies that are discarded
	// after scoring, so atom state persistence doesn't matter.
	// The serial output path (cluster.cpp:240→302, top.cpp:1242→1287) calls
	// ic2cf precisely to populate atoms[] with the final docked Cartesian pose,
	// then immediately passes atoms[] to write_pdb.  Restoring atoms on normal
	// exit would silently write the pre-call (initial) structure — every cluster
	// output PDB would be a copy of the start conformation, not the docked pose.
	// Atom restore belongs only on the OOB/penalty early-exit path above, where
	// we must leave a clean baseline for the next chromosome evaluation.
	FA->ori[0] = ori_save[0];
	FA->ori[1] = ori_save[1];
	FA->ori[2] = ori_save[2];

	return cf;

}

#ifdef _WIN32
double get_apparent_cf_evalue(cfstr* cf) {
#else
	double get_apparent_cf_evalue(cfstr* cf) {
#endif
		return cf->com + cf->wal + cf->sas + cf->elec + cf->hbond + cf->gist_desolv + cf->metal_coord;
	}

#ifdef _WIN32
	double get_cf_evalue(cfstr* cf) {
#else
		double get_cf_evalue(cfstr* cf) {
#endif
			return cf->com + cf->wal + cf->sas + cf->con + cf->elec + cf->hbond + cf->gist_desolv + cf->metal_coord;
		}
