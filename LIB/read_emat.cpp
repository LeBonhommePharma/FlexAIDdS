#include <sstream>
#include <string>
#include <vector>
#include <cstring>  // memcpy

#include "flexaid.h"
#include "fileio.h"

void read_emat(FA_Global* FA, char* emat_file)
{
    FILE* infile_ptr;
    char buffer[100];

    infile_ptr=NULL;
    if (!OpenFile_B(emat_file,"r",&infile_ptr)){
        fprintf(stderr,"ERROR: Could not read input file: %s\n", emat_file);
        Terminate(8);
    }

    std::vector<std::string> lines;
    std::string pairwiseline = "";

    // builds a string from a buffered line
    while(fgets(buffer,sizeof(buffer),infile_ptr) != NULL){
        pairwiseline += std::string(buffer);
        if(pairwiseline.find("\n") != std::string::npos){
            lines.push_back(pairwiseline);
            pairwiseline = "";
        }
    }

    printf("read %d lines in <%s>\n", (int)lines.size(), emat_file);

    float z = zero(1.0/2.0f, 1.0/2.0f, (float)(-(int)lines.size()));
    if(fabs(z - (float)((int)z)) > 0.001) {
        fprintf(stderr,"ERROR: Number of lines read in energy matrix file <%s> is incorrect (%d)\n",
                emat_file, (int)lines.size());
        Terminate(12);
    }

    FA->ntypes = (int)(z + 0.001);
    printf("number of atom types: %d\n", FA->ntypes);

    FA->energy_matrix = (struct energy_matrix*)malloc(FA->ntypes*FA->ntypes*sizeof(struct energy_matrix));
    if(!FA->energy_matrix){
        fprintf(stderr,"ERROR: could not allocate memory for energy_matrix\n");
        Terminate(2);
    }

    for(int i=0;i<FA->ntypes;i++){
        for(int j=i;j<FA->ntypes;j++){
            std::string ori_line, line = *lines.begin();
            ori_line = line;
            lines.erase(lines.begin());

            if(line.find("=") != std::string::npos){
                line = line.substr(line.find("=")+1);
            }

            // Start of replacement
            const char* whitespace = " \t\n";

            // Trim leading whitespace
            size_t first = line.find_first_not_of(whitespace);
            if (std::string::npos == first) {
                line = "";
            } else {
                // Trim trailing whitespace
                size_t last = line.find_last_not_of(whitespace);
                line = line.substr(first, (last - first + 1));
            }

            // Split string into tokens
            std::vector<std::string> values;
            std::stringstream ss(line);
            std::string token;
            while (ss >> token) {
                values.push_back(token);
            }
            // End of replacement

            FA->energy_matrix[i*FA->ntypes+j].type1 = i+1;
            FA->energy_matrix[i*FA->ntypes+j].type2 = j+1;
            FA->energy_matrix[j*FA->ntypes+i].type1 = j+1;
            FA->energy_matrix[j*FA->ntypes+i].type2 = i+1;
            // A2: zero-init flat arrays for all cases; populated below for weight=0
            FA->energy_matrix[i*FA->ntypes+j].flat_n = 0;
            FA->energy_matrix[i*FA->ntypes+j].flat_x = NULL;
            FA->energy_matrix[i*FA->ntypes+j].flat_y = NULL;
            FA->energy_matrix[i*FA->ntypes+j].flat_slope = NULL;
            FA->energy_matrix[j*FA->ntypes+i].flat_n = 0;
            FA->energy_matrix[j*FA->ntypes+i].flat_x = NULL;
            FA->energy_matrix[j*FA->ntypes+i].flat_y = NULL;
            FA->energy_matrix[j*FA->ntypes+i].flat_slope = NULL;

            if(values.size() == 1){
                FA->energy_matrix[i*FA->ntypes+j].weight = 1;
                FA->energy_matrix[j*FA->ntypes+i].weight = 1;

                struct energy_values* weightval = (struct energy_values*)malloc(sizeof(struct energy_values));
                if(!weightval){
                    fprintf(stderr,"ERROR: could not allocate memory for weightval\n");
                    Terminate(2);
                }

                weightval->x = -1;
                weightval->y = atof((*values.begin()).c_str());
                weightval->next_value = NULL;

                FA->energy_matrix[i*FA->ntypes+j].energy_values = weightval;
                FA->energy_matrix[j*FA->ntypes+i].energy_values = weightval;

            }else if(values.size() % 2 == 0){
                FA->energy_matrix[i*FA->ntypes+j].weight = 0;
                FA->energy_matrix[j*FA->ntypes+i].weight = 0;

                struct energy_values* xyval_prev = NULL;

                for(std::vector<std::string>::iterator it=values.begin(); it!=values.end(); it+=2){
                    std::vector<std::string>::iterator xit = it;
                    std::vector<std::string>::iterator yit = it+1;

                    struct energy_values* xyval = (struct energy_values*)malloc(sizeof(struct energy_values));
                    if(!xyval){
                        fprintf(stderr,"ERROR: could not allocate memory for xyval\n");
                        Terminate(2);
                    }

                    xyval->x = atof((*xit).c_str());
                    xyval->y = atof((*yit).c_str());
                    xyval->next_value = NULL;

                    // multiply by a factor otherwise only conformer-driven
                    //xyval->y *= 10.0;

                    /*
                      printf("new entry[%d][%d]: x=%.3f y=%.3f\n", i+1, j+1,
                             xyval->x, xyval->y);
                    */

                    // second or more xy values
                    if(xyval_prev != NULL)
                        xyval_prev->next_value = xyval;
                    else {
                        FA->energy_matrix[i*FA->ntypes+j].energy_values = xyval;
                        FA->energy_matrix[j*FA->ntypes+i].energy_values = xyval;
                    }

                    xyval_prev = xyval;
                }

                // A2: build flat arrays for fast get_yval (no linked-list walk per eval)
                {
                    int fn = (int)(values.size() / 2);
                    // layout: [flat_x | flat_y | flat_slope], one contiguous malloc
                    float* buf_ij = (float*)malloc(3 * fn * sizeof(float));
                    if(buf_ij) {
                        float* fx = buf_ij;
                        float* fy = buf_ij + fn;
                        float* fs = buf_ij + 2*fn;
                        int k = 0;
                        for(struct energy_values* ev = FA->energy_matrix[i*FA->ntypes+j].energy_values;
                            ev; ev = ev->next_value, ++k) {
                            fx[k] = ev->x;  fy[k] = ev->y;
                        }
                        for(int s = 0; s < fn-1; ++s)
                            fs[s] = (fy[s+1] - fy[s]) / (fx[s+1] - fx[s]);
                        fs[fn-1] = 0.0f;  // unused (last segment has no right neighbor)
                        FA->energy_matrix[i*FA->ntypes+j].flat_n = fn;
                        FA->energy_matrix[i*FA->ntypes+j].flat_x = fx;
                        FA->energy_matrix[i*FA->ntypes+j].flat_y = fy;
                        FA->energy_matrix[i*FA->ntypes+j].flat_slope = fs;
                        // symmetric entry gets its own copy
                        float* buf_ji = (float*)malloc(3 * fn * sizeof(float));
                        if(buf_ji) {
                            memcpy(buf_ji, buf_ij, 3 * fn * sizeof(float));
                            FA->energy_matrix[j*FA->ntypes+i].flat_n = fn;
                            FA->energy_matrix[j*FA->ntypes+i].flat_x = buf_ji;
                            FA->energy_matrix[j*FA->ntypes+i].flat_y = buf_ji + fn;
                            FA->energy_matrix[j*FA->ntypes+i].flat_slope = buf_ji + 2*fn;
                        }
                    }
                }
            }else{
                fprintf(stderr,"ERROR: invalid number of xy-values for atom pairwise %d-%d\n", i, j);
                Terminate(12);
            }
        }
    }

    CloseFile_B(&infile_ptr,"r");
}