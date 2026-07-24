#include "flexaid.h"
#include "fileio.h"
#include "maps.hpp"
#include <cstdlib>   // getenv, atof

/*******************************************/
/*  generates a grid from the spheres of
    loccen or spheres combination from
    locclf (using a cleft)                 */
/*******************************************/

gridpoint* generate_grid(FA_Global* FA,sphere* spheres, atom* atoms, resid* residue){

	float sqrrad;
	float c[3], min[3], max[3];
	gridpoint* cleftgrid = NULL;
	std::map<GridKey,int> cleftgrid_map;

	cleftgrid = (gridpoint*)malloc(FA->MIN_CLEFTGRID_POINTS*sizeof(gridpoint));
	if (cleftgrid == NULL){
		fprintf(stderr,"ERROR: memory allocation error for cleftgrid\n");
		Terminate(2);
	}

    cleftgrid[0].number = 0;
    for (int j=0;j<3;j++) cleftgrid[0].coor[j] = atoms[residue[FA->res_cnt].gpa[0]].coor[j];
    cleftgrid[0].dis = atoms[residue[FA->res_cnt].gpa[0]].dis;
    cleftgrid[0].ang = atoms[residue[FA->res_cnt].gpa[0]].ang;
    cleftgrid[0].dih = atoms[residue[FA->res_cnt].gpa[0]].dih;

	// ── P3 finer-grid lever: FLEXAIDDS_GRID_SPACING (Å) ───────────────────
	// Env-gated override of the grid spacer length used to tile the cleft.
	// Default-preserving: unset (or <= 0) leaves FA->spacer_length untouched,
	// so the historical 0.375 Å default (config optimization.grid_spacing) is
	// byte-stable. A finer spacing resolves sub-0.1 Å IC-only native basins
	// (e.g. 1K3U's ~0.078 Å basin) that a coarse grid cannot represent.
	// Mutating the shared FA->spacer_length keeps slice_grid / partition_grid
	// (which read the same field) consistent with the finer grid. GridKey snaps
	// to milliangstrom, so spacings down to 0.05 Å dedup correctly.
	// This is a SEARCH-COVERAGE lever: expected to pay off only combined with
	// anti-collapse (P1). On its own it just enlarges the grid (more memory /
	// vertices), so it is clamped to [0.05, 2.0] Å to bound the explosion.
	{
		static const char* gs_env = std::getenv("FLEXAIDDS_GRID_SPACING");
		if (gs_env != NULL && gs_env[0] != '\0'){
			float gs = (float)std::atof(gs_env);
			if (gs > 0.0f){
				if (gs < 0.05f) gs = 0.05f;
				if (gs > 2.0f)  gs = 2.0f;
				printf("[GRID-SPACING] override FA->spacer_length %.3f -> %.3f "
				       "(FLEXAIDDS_GRID_SPACING)\n", FA->spacer_length, gs);
				FA->spacer_length = gs;
			}
		}
	}

	printf("will build a grid with spacing %.3f\n", FA->spacer_length);

	FA->num_grd = 1; // set counter to 1 because 0 is the INI conformation of the ligand
	while(spheres != NULL){

		if ( (float)( 1.0 / FA->spacer_length ) - (float)((int)( 1.0 / FA->spacer_length )) > 0.001 ){
			min[0] = (float)((int)( (spheres->center[0] - spheres->radius) / FA->spacer_length )) * FA->spacer_length;
			min[1] = (float)((int)( (spheres->center[1] - spheres->radius) / FA->spacer_length )) * FA->spacer_length;
			min[2] = (float)((int)( (spheres->center[2] - spheres->radius) / FA->spacer_length )) * FA->spacer_length;
			max[0] = (float)((int)( (spheres->center[0] + spheres->radius) / FA->spacer_length ) + 1.0) * FA->spacer_length;
			max[1] = (float)((int)( (spheres->center[1] + spheres->radius) / FA->spacer_length ) + 1.0) * FA->spacer_length;
			max[2] = (float)((int)( (spheres->center[2] + spheres->radius) / FA->spacer_length ) + 1.0) * FA->spacer_length;
		}else{
			min[0] = (float)((int)( (spheres->center[0] - spheres->radius - FA->spacer_length )));
			min[1] = (float)((int)( (spheres->center[1] - spheres->radius - FA->spacer_length )));
			min[2] = (float)((int)( (spheres->center[2] - spheres->radius - FA->spacer_length )));
			max[0] = (float)((int)( (spheres->center[0] + spheres->radius + FA->spacer_length ) + 1.0));
			max[1] = (float)((int)( (spheres->center[1] + spheres->radius + FA->spacer_length ) + 1.0));
			max[2] = (float)((int)( (spheres->center[2] + spheres->radius + FA->spacer_length ) + 1.0));
		}

		c[0] = min[0];
		c[1] = min[1];
		c[2] = min[2];

		sqrrad = spheres->radius * spheres->radius;

		while(c[2] < max[2]){
			while(c[1] < max[1]){
				while(c[0] < max[0]){

					if(sqrdist(spheres->center,c) < sqrrad){

						GridKey key(c);

						if(cleftgrid_map.find(key) == cleftgrid_map.end()){
							if (FA->num_grd==FA->MIN_CLEFTGRID_POINTS){
								FA->MIN_CLEFTGRID_POINTS *= 2;

								cleftgrid = (gridpoint*)realloc(cleftgrid,FA->MIN_CLEFTGRID_POINTS*sizeof(gridpoint));
								if (cleftgrid == NULL){
									fprintf(stderr,"ERROR: memory reallocation error for cleftgrid\n");
									Terminate(2);
								}
							}

							memset(&cleftgrid[FA->num_grd], 0, sizeof(gridpoint));
							cleftgrid[FA->num_grd].coor[0] = c[0];
							cleftgrid[FA->num_grd].coor[1] = c[1];
							cleftgrid[FA->num_grd].coor[2] = c[2];

							FA->num_grd++;
							cleftgrid_map.insert(std::pair<GridKey,int>(key, FA->num_grd));

						}

					}

					c[0] += FA->spacer_length;
				}

				c[0] = min[0];
				c[1] += FA->spacer_length;
			}

			c[0] = min[0];
			c[1] = min[1];
			c[2] += FA->spacer_length;

		}

		spheres = spheres->prev;
	}

	printf("built a grid with %d vertices\n", FA->num_grd - 1);

	return cleftgrid;
}
