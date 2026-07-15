#include "flexaid.h"
#include "fileio.h"
#include "ensemble_pipeline.h"

#include <cstdlib>
#include <cstring>

sphere* read_spheres(char filename[]){

	FILE* file_ptr;             // read file stream
	char buffer[128];           // buffer to read line (was 81; need ≥66 for B-factor radius)
	char field[7];              // 6 first character of the line
	char coor[9];               // coordinate
	char radius[8];             // radius field
	int i,j,k;                  // dummy counters
	sphere* spheres = NULL;     // spheres list
	int n_rejected = 0;
	
	file_ptr=NULL;
	if (!OpenFile_B(filename,"r",&file_ptr))
		Terminate(8);
  
	while (fgets(buffer, sizeof(buffer), file_ptr) != NULL){
		
		//0         1         2         3         4         5         6         7
		//0123456789012345678901234567890123456789012345678901234567890123456789
		//ATOM   1102  C   SPH Z   1      34.069  28.877   7.194  1.00  1.51 
		for (i=0;i<6;i++) field[i] = buffer[i];
		field[6] = '\0';
    
		if (strncmp(buffer, "ATOM  ", 6) == 0 || strncmp(buffer, "HETATM", 6) == 0){
			// Reproducibility: reject short/malformed lines (silent r=0 → empty grid).
			const size_t len = std::strlen(buffer);
			if (len < 66) {
				++n_rejected;
				continue;
			}

			sphere* _sphere;
			_sphere = (sphere*)malloc(sizeof(sphere));
			if(_sphere == NULL){
				fprintf(stderr,"ERROR: memory allocation error for spheres (LOCCLF)\n");
				Terminate(2);
			}
			
			
			for (j=0;j<=2;j++){
				k=0;
				for (i=30+j*8;i<30+(j+1)*8;i++)
					coor[k++] = buffer[i];
				coor[8] = '\0';
				
				if (sscanf(coor, "%f", &_sphere->center[j]) != 1) {
					free(_sphere);
					++n_rejected;
					_sphere = NULL;
					break;
				}
			}
			if (!_sphere) continue;
			
			// B-factor columns 61-66 hold radius (GetCleft / CleftDetector convention).
			for(i=0; i<5; i++)
				radius[i] = buffer[i+61];
			radius[5] = '\0';

			char* endp = nullptr;
			const float r = static_cast<float>(std::strtod(radius, &endp));
			if (endp == radius || !ensemble::valid_sphere_radius(r)) {
				free(_sphere);
				++n_rejected;
				continue;
			}
			_sphere->radius = r;

			_sphere->prev = spheres;
			spheres = _sphere;
		}
	}

	if (n_rejected > 0) {
		fprintf(stderr,
			"[read_spheres] rejected %d malformed/zero-radius sphere line(s) in %s\n",
			n_rejected, filename);
	}
	if (spheres == NULL) {
		fprintf(stderr,
			"[read_spheres] ERROR: no valid spheres in %s (empty grid would not be reproducible)\n",
			filename);
	}
  
	return spheres;
}
