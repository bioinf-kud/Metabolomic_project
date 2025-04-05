#include <stdio.h>
#include <stdlib.h>
#include <string.h>
struct reactPair {
    char reactant[60];
    char product[60];
};
int main() {
    char names_dir[]="/Users/sunkai/Desktop/metabolomics_project/c_project/data/find_reaction/basic_data/finding_cofactors/metabolite_names.tsv";
    char unpair_dir[]="/Users/sunkai/Desktop/metabolomics_project/c_project/data/find_reaction/basic_data/finding_cofactors/unpaired.tsv";
    char start_dir[]="/Users/sunkai/Desktop/metabolomics_project/c_project/data/find_reaction/basic_data/finding_cofactors/start.tsv";
    char end_dir[]="/Users/sunkai/Desktop/metabolomics_project/c_project/data/find_reaction/basic_data/finding_cofactors/end.tsv";
    char outunpair_dir[]="/Users/sunkai/Desktop/metabolomics_project/c_project/data/find_reaction/basic_data/finding_cofactors/unpaired_id.tsv";
    char outpair_dir[]="/Users/sunkai/Desktop/metabolomics_project/c_project/data/find_reaction/basic_data/finding_cofactors/paired_id.tsv";
    FILE *names_file = fopen(names_dir, "r");
    FILE *unpair_file = fopen(unpair_dir, "r");
    FILE *start_file = fopen(start_dir, "r");
    FILE *end_file = fopen(end_dir, "r");
    FILE *outunpair_file = fopen(outunpair_dir, "w+");
    FILE *outpair_file = fopen(outpair_dir, "w+");
    char metabolite[8399][60];
    char metabolite_simple[8399][60];
    char metabolite_compart[8399];
    struct reactPair reactpair[52];
    char unpair[11][60];
    for (int i = 0; i < 11; i++) {
       fscanf(unpair_file, "%s", unpair[i]);
    }
    for (int i = 0; i < 8399; i++) {
        fscanf(names_file, "%s", metabolite[i]);
    }
    for (int i = 0; i < 52; i++) {
        fscanf(start_file, "%s", reactpair[i].reactant);
        fscanf(end_file, "%s", reactpair[i].product);
    }
    for (int i = 0; i < 8399; i++) {
        char *start_bracket = strchr(metabolite[i], '[');
        if (start_bracket != NULL) {
            int prefix_length = start_bracket - metabolite[i];
            strncpy(metabolite_simple[i], metabolite[i], prefix_length);
            metabolite_simple[i][prefix_length] = '\0'; 
            metabolite_compart[i] = *(start_bracket + 1); 
        } else {
            strcpy(metabolite_simple[i], metabolite[i]);
            metabolite_compart[i] = '\0'; 
        }
    }
    for(int i=0; i<11; i++){
        for(int j=0; j<8399; j++){
            if(strcmp(unpair[i], metabolite_simple[j])==0){
                fprintf(outunpair_file, "%d\t%s\n", j, metabolite[j]);
            }
        }
    }
    for(int i=0; i<52; i++){
        char tempcomp;
        for(int j=0; j<8399; j++)
        {
            if(strcmp(reactpair[i].reactant, metabolite_simple[j])==0){
                tempcomp = metabolite_compart[j];
                for(int k=0; k<8399; k++){
                    if(strcmp(reactpair[i].product, metabolite_simple[k])==0 && tempcomp == metabolite_compart[k]){
                        fprintf(outpair_file, "%d\t%s\t%d\t%s\n", j, metabolite[j], k, metabolite[k]);
                    }
                }
            }
        }
    }
    fclose(names_file);
    fclose(unpair_file);
    fclose(start_file);
    fclose(end_file);
    fclose(outunpair_file);
    fclose(outpair_file);
    printf("Finished!\n");
    return 0;
}
