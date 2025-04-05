#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define REACTION_COUNT 13543 
#define METABOLITE_COUNT 8399
#define REACTION_TABLE_SIZE 13553 // A prime number greater than REACTION_COUNT  
#define METABOLITE_TABLE_SIZE 8419  // A prime number greater than METABOLITE_COUNT
#define NAME_SIZE 60
typedef struct metabolite_adjacency_list {
    int adjacent;
    int reaction_number;
    int*reactions;
    int*reactions_direction;
    double enzyme_effeciency;
} metabolite_adjacency_list;
// Hash table node
typedef struct Node {
    char *key;
    int value;
    struct Node *next;
} Node;
// Dynamic array
typedef struct {
    double *rxn_data;
    int *id;      
    int size;       
    int capacity;   
} DynamicArray;
typedef struct {
    int sub_id;      
    int prod_id; 
} cof_pair;
// Hash tables
Node *reaction_hash_table[REACTION_TABLE_SIZE];
Node *metabolite_hash_table[METABOLITE_TABLE_SIZE];
metabolite_adjacency_list adjacency_node[METABOLITE_COUNT][METABOLITE_COUNT];

unsigned int hash(const char *key, int table_size);// Hash function
void insert(Node *hash_table[], const char *key, int value, int table_size);// Insert a key-value pair into the hash table
int lookup(Node *hash_table[], const char *key, int table_size);// Lookup a value by key in the hash table
void free_table(Node *hash_table[], int table_size);// Free the memory used by the hash table
int setup_reaction_hash_table(Node *hash_table[],char* file_path,char** reaction_names);
int setup_metabolite_hash_table(Node *hash_table[],char* file_path,char** metabolite_names);// Setup the reaction hash table
DynamicArray* create_array(int initial_capacity);// initialize a dynamic array 
void add_element(DynamicArray *arr, int id,double value);// add an element to the dynamic array
void remove_element(DynamicArray *arr, int index);// remove an element from the dynamic array by index
void free_array(DynamicArray *arr);// free space allocated for the dynamic array
int read_cofactors(char* unp_file_path,char* p_file_path,int *unpaired_id, cof_pair *paired_id);//read the cofactors from the file
int scale_reaction(char* rxn_met_file_path,char *rxn_with_data_file_path,char* out_path,char* out_num_path,int *unpaired_id, cof_pair *paired_id, int unp_num, int p_num, char **reaction_names, char **metabolite_names);//scale the reaction
int delete_cofactor(int *unpaired_id, cof_pair *paired_id, int unp_num, int p_num, DynamicArray *sub_arr, DynamicArray *prod_arr);//delete the cofactors
int delete_coa(DynamicArray *sub_arr, DynamicArray *prod_arr, char **metabolite_names);//delete the coa related cofactors
int manual_delete_output(int reactnum,FILE*out,FILE* outnum,DynamicArray *sub_arr, DynamicArray *prod_arr, char **metabolite_names,char **reaction_names);//delete the coa related cofactors
int output_reaction(int reactnum,FILE*out,FILE* outnum, char **reaction_names, char **metabolite_names,DynamicArray *sub_arr, DynamicArray *prod_arr);//output the reaction
void remove_last_three_chars(char *name);//remove the last three characters from the string
int ends_with_coa(const char *str);//check if the string ends with "coa"
int main() {
    // Reaction names array (ID -> name)
    const char *reaction_names[REACTION_COUNT];
    const char *metabolite_names[METABOLITE_COUNT];
    const int unp_num=76;
    const int p_num=185;
    int unpaired_id[unp_num];
    cof_pair paired_id[p_num];
    char rxn_list_dir[] = "../data/scale_reaction/input/reaction_and_metabolite_info/rxns_list.tsv";
    char met_list_dir[] = "../data/scale_reaction/input/reaction_and_metabolite_info/metabolites_list.tsv";
    char rxn_met_dir[] = "../data/scale_reaction/input/reaction_and_metabolite_info/reaction_mat.tsv";
    char rxn_with_data_dir[] = "../data/scale_reaction/input/target_reaction_id/reaction_with_data_num.tsv";
    char out_num_dir[] = "../data/scale_reaction/output/reaction_substrate_num_list.tsv";
    char out_dir[] = "../data/scale_reaction/output/reaction_dataframe.txt";
    char unpair_cof_dir[] ="../data/scale_reaction/input/cofactor_data/unpaired_id.tsv";
    char pair_cof_dir[] ="../data/scale_reaction/input/cofactor_data/paired_id.tsv";
    // Setup hash tables for reactions and metabolites
    printf("Setting up reaction hash table...\n");
    int a=setup_reaction_hash_table(reaction_hash_table,rxn_list_dir,(char**)reaction_names);
    if (a == -1) {
        fprintf(stderr, "Error setting up reaction hash table\n");
        return 1;
    }
    printf("Reaction hash table setup complete.\n");
    printf("Setting up metabolite hash table...\n");
    int b=setup_metabolite_hash_table(metabolite_hash_table,met_list_dir,(char**)metabolite_names);
    if (b == -1) {
        fprintf(stderr, "Error setting up metabolite hash table\n");
        return 1;
    }
    printf("Metabolite hash table setup complete.\n");
    // Scale reactions
    printf("Reading cofactors...\n");
    read_cofactors(unpair_cof_dir,pair_cof_dir,unpaired_id,paired_id);
    printf("Scaling reactions...\n");
    scale_reaction(rxn_met_dir,rxn_with_data_dir,out_dir,out_num_dir,unpaired_id,paired_id,unp_num,p_num,(char**)reaction_names,(char**)metabolite_names);
    // Free memory
    free_table(reaction_hash_table, REACTION_TABLE_SIZE);
    free_table(metabolite_hash_table, METABOLITE_TABLE_SIZE);
    for (int i = 0; i < REACTION_COUNT; i++) {
        free((void *)reaction_names[i]);
    }
    for (int i = 0; i < METABOLITE_COUNT; i++) {
        free((void *)metabolite_names[i]);
    }
    return 0;
}

unsigned int hash(const char *key, int table_size) {
    unsigned int hash = 0;
    while (*key) {
        hash = (hash << 5) + *key++;
    }
    return hash % table_size;
}
void insert(Node *hash_table[], const char *key, int value, int table_size) {
    unsigned int index = hash(key, table_size);
    Node *new_node = malloc(sizeof(Node));
    new_node->key = strdup(key);
    new_node->value = value;
    new_node->next = hash_table[index];
    hash_table[index] = new_node;
}
int lookup(Node *hash_table[], const char *key, int table_size) {
    unsigned int index = hash(key, table_size);
    Node *current = hash_table[index];
    while (current) {
        if (strcmp(current->key, key) == 0) {
            return current->value;
        }
        current = current->next;
    }
    return -1; // Not found
}
void free_table(Node *hash_table[], int table_size) {
    for (int i = 0; i < table_size; i++) {
        Node *current = hash_table[i];
        while (current) {
            Node *temp = current;
            current = current->next;
            free(temp->key);
            free(temp);
        }
    }
}
int setup_reaction_hash_table(Node *hash_table[],char* file_path,char** reaction_names) {
    FILE *file = fopen(file_path, "r");
    if (!file) {
        perror("Cannot open reaction file");
        return -1;
    }

    for(int i = 0; i < REACTION_COUNT; i++) {
        char line[NAME_SIZE];
        if (fgets(line, sizeof(line), file) == NULL) {
            break;
        }
        char *name = malloc(NAME_SIZE);
        reaction_names[i] = name;
        sscanf(line, "%s", name);
        insert(hash_table, name, i, REACTION_TABLE_SIZE);
    }
    fclose(file);
    return 0;
}
int setup_metabolite_hash_table(Node *hash_table[],char* file_path,char** metabolite_names) {
    FILE *file = fopen(file_path, "r");
    if (!file) {
        perror("Cannot open metabolite file");
        return -1;
    }

    for(int i = 0; i < METABOLITE_COUNT; i++) {
        char line[NAME_SIZE];
        if (fgets(line, sizeof(line), file) == NULL) {
            break;
        }
        char *name = malloc(NAME_SIZE);
        metabolite_names[i] = name;
        sscanf(line, "%s", name);
        insert(hash_table, name, i, METABOLITE_TABLE_SIZE);
    }
    fclose(file);
    return 0;
}
DynamicArray* create_array(int initial_capacity) {
    DynamicArray *arr = malloc(sizeof(DynamicArray));
    arr->id = malloc(sizeof(int) * initial_capacity);
    arr->rxn_data = malloc(sizeof(double) * initial_capacity);
    arr->size = 0;
    arr->capacity = initial_capacity;
    return arr;
}
void add_element(DynamicArray *arr, int id,double value) {
    if (arr->size == arr->capacity) {
        // expand the array
        arr->capacity *= 2;
        arr->id = realloc(arr->id, sizeof(int) * arr->capacity);
        arr->rxn_data = realloc(arr->rxn_data, sizeof(double) * arr->capacity);
    }
    arr->id[arr->size] = id;
    arr->rxn_data[arr->size] = value;
    arr->size++;
}
void remove_element(DynamicArray *arr, int index) {
    if (index < 0 || index >= arr->size) {
        printf("Index out of bounds\n");
        return;
    }
    for (int i = index; i < arr->size - 1; i++) {
        arr->id[i] = arr->id[i + 1];
        arr->rxn_data[i] = arr->rxn_data[i + 1];
    }
    arr->size--;
}
void free_array(DynamicArray *arr) {
    free(arr->id);
    free(arr->rxn_data);
    free(arr);
}
int read_cofactors(char* unp_file_path,char* p_file_path,int *unpaired_id, cof_pair *paired_id){
    printf("Reading unpaired cofactors...\n");
    FILE *f = fopen(unp_file_path, "r");
    if (!f) {
        perror("Cannot open unpaired cofactors file");
        return -1;
    }
    for (int i = 0; i < 76; i++) {
        fscanf(f, "%d", &unpaired_id[i]);
    }
    fclose(f);
    printf("Reading paired cofactors...\n");
    f = fopen(p_file_path, "r");
    if (!f) {
        perror("Cannot open paired cofactors file");
        return -1;
    }
    for (int i = 0; i < 185; i++) {
        fscanf(f, "%d", &paired_id[i].sub_id);
        fscanf(f, "%d", &paired_id[i].prod_id);
    }
    fclose(f);
    return 0;
}
int scale_reaction(char* rxn_met_file_path,char *rxn_with_data_file_path,char* out_path,char* out_num_path,int *unpaired_id, cof_pair *paired_id, int unp_num, int p_num, char **reaction_names, char **metabolite_names){
    FILE *f = fopen(rxn_met_file_path, "r");
    FILE *f2 = fopen(rxn_with_data_file_path, "r");
    FILE *out = fopen(out_path, "w+");
    FILE *out_num = fopen(out_num_path, "w+");
    if (!f || !f2 || !out|| !out_num) {
        perror("Cannot open file");
        return -1;
    }
    // Read and scale the reaction data
    printf("Reading and scaling reactions...\n");
    double metabolites[METABOLITE_COUNT];
    int cur_reaction = 0;
    int rcnt = 0;
    fscanf(f2, "%d", &cur_reaction);
    for(int i = 0; i < REACTION_COUNT; i++) {
        for(int j = 0; j < METABOLITE_COUNT; j++) {
            fscanf(f, "%lf", &metabolites[j]);
        }
        if(i == cur_reaction) {
            rcnt++;
            printf("Processing reaction #%d, %s...\n",rcnt,reaction_names[cur_reaction]);
            int subcnt = 0;
            int prodcnt = 0;
            for(int j = 0; j < METABOLITE_COUNT; j++) {
                if(metabolites[j] > 0) {
                    subcnt++;
                } else if(metabolites[j] < 0) {
                    prodcnt++;
                }
            }
            DynamicArray *sub_array = create_array(subcnt);
            DynamicArray *prod_array = create_array(prodcnt);
            for(int j = 0; j < METABOLITE_COUNT; j++) {
                if(metabolites[j] > 0) {
                    add_element(sub_array, j, metabolites[j]);
                } else if(metabolites[j] < 0) {
                    add_element(prod_array, j, metabolites[j]);
                }
            }
            if(sub_array->size > 1 || prod_array->size > 1)
                delete_cofactor(unpaired_id, paired_id, unp_num, p_num, sub_array, prod_array);
            if(sub_array->size > 1 || prod_array->size > 1) 
                delete_coa(sub_array, prod_array, metabolite_names);        
            if(sub_array->size > 1 || prod_array->size > 1) {
                manual_delete_output(cur_reaction,out, out_num,sub_array, prod_array, metabolite_names, reaction_names);
                output_reaction(cur_reaction, out, out_num, reaction_names, metabolite_names, sub_array, prod_array);
            }
            else{
                output_reaction(cur_reaction, out, out_num, reaction_names, metabolite_names, sub_array, prod_array);
            }
        
            free_array(sub_array);
            free_array(prod_array);
            if (fscanf(f2, "%d", &cur_reaction) == EOF) {
                break; 
            }
        }
    }
    fclose(out);
    fclose(f);
    fclose(f2);
    return 0;
}
int delete_cofactor(int *unpaired_id, cof_pair *paired_id, int unp_num, int p_num, DynamicArray *sub_arr, DynamicArray *prod_arr) {
    for (int i = 0; i < unp_num; i++) {
        for (int j = 0; j < sub_arr->size; j++) {
            if (sub_arr->id[j] == unpaired_id[i]) {
                remove_element(sub_arr, j);
                break;
            }
        }
    }
    for (int i = 0; i < unp_num; i++) {
        for (int j = 0; j < prod_arr->size; j++) {
            if (prod_arr->id[j] == unpaired_id[i]) {
                remove_element(prod_arr, j);
                break;
            }
        }
    }
    if(sub_arr->size == 1 || prod_arr->size == 1){
        return 0; 
    }
    for (int i = 0; i < p_num; i++) {
        for (int j = 0; j < sub_arr->size; j++) {
            if (sub_arr->id[j] == paired_id[i].sub_id) {
                for(int k = 0; k < prod_arr->size; k++) {
                    if (prod_arr->id[k] == paired_id[i].prod_id) {
                        remove_element(sub_arr, j);
                        remove_element(prod_arr, k);
                        if(sub_arr->size == 1|| prod_arr->size == 1){
                            return 0; 
                        }
                        break;
                    }
                }
            }
        }
        for (int j = 0; j < prod_arr->size; j++) {
            if (prod_arr->id[j] == paired_id[i].sub_id) {
                for(int k = 0; k < sub_arr->size; k++) {
                    if (sub_arr->id[k] == paired_id[i].prod_id) {
                        remove_element(prod_arr, j);
                        remove_element(sub_arr, k);
                        if(sub_arr->size == 1|| prod_arr->size == 1){
                            return 0; 
                        }
                        break;
                    }
                }
            }
        }  
    }
    return 0;
}
int delete_coa(DynamicArray *sub_arr, DynamicArray *prod_arr, char**metabolite_names) {
    if(sub_arr->size == 1){// R-transfer reaction
        if(prod_arr->size > 1){
            for(int i = 0; i < prod_arr->size; i++){
                char name[NAME_SIZE];
                strcpy(name, metabolite_names[prod_arr->id[i]]);
                remove_last_three_chars(name);
                if(strcmp(name, "coa") == 0){
                    remove_element(prod_arr, i);
                    return 0;
            }
        }
    }
    } else if(prod_arr->size == 1){// Coa-linking reaction
        if(sub_arr->size > 1){
            for(int i = 0; i < sub_arr->size; i++){
                char name[NAME_SIZE];
                strcpy(name, metabolite_names[sub_arr->id[i]]);
                remove_last_three_chars(name);
                if(strcmp(name, "coa") == 0){
                    remove_element(sub_arr, i);
                    return 0;
                }
            }
        }
    }
    int subcnt = 0;
    int prodcnt = 0;
    for(int i = 0; i < sub_arr->size; i++){
        char sub_name[NAME_SIZE];
        strcpy(sub_name, metabolite_names[sub_arr->id[i]]);
        remove_last_three_chars(sub_name);
        if(ends_with_coa(sub_name) == 1)
            subcnt++;
    }
    for(int i = 0; i < prod_arr->size; i++){
        char prod_name[NAME_SIZE];
        strcpy(prod_name, metabolite_names[prod_arr->id[i]]);
        remove_last_three_chars(prod_name);
        if(ends_with_coa(prod_name) == 1)
            prodcnt++;
    }
    if(subcnt == 1 && prodcnt == 1){//remove paired coa
        for(int i = 0; i < sub_arr->size; i++){
            char sub_name[NAME_SIZE];
            strcpy(sub_name, metabolite_names[sub_arr->id[i]]);
            remove_last_three_chars(sub_name);
            if(ends_with_coa(sub_name) == 1){
                for(int j = 0; j < prod_arr->size; j++){
                    char prod_name[NAME_SIZE];
                    strcpy(prod_name, metabolite_names[prod_arr->id[j]]);
                    remove_last_three_chars(prod_name);
                    if(ends_with_coa(prod_name) == 1){
                        remove_element(sub_arr, i);
                        remove_element(prod_arr, j);
                        return 0;
                    }
                }
            }
        }
    }
    if(sub_arr->size >1 && prod_arr->size > 1){//remove accoa-coa and malcoa-coa
        for(int i = 0; i < sub_arr->size; i++){
            char sub_name[NAME_SIZE];
            strcpy(sub_name, metabolite_names[sub_arr->id[i]]);
            remove_last_three_chars(sub_name);
            if(strstr(sub_name, "coa") != NULL){
                for(int j = 0; j < prod_arr->size; j++){
                    char prod_name[NAME_SIZE];
                    strcpy(prod_name, metabolite_names[prod_arr->id[j]]);
                    remove_last_three_chars(prod_name);
                    if(strstr(prod_name, "accoa") != NULL){
                        remove_element(sub_arr, i);
                        remove_element(prod_arr, j);
                        return 0;
                    }
                }
            }
        }
        for(int i = 0; i < prod_arr->size; i++){
            char prod_name[NAME_SIZE];
            strcpy(prod_name, metabolite_names[prod_arr->id[i]]);
            remove_last_three_chars(prod_name);
            if(strstr(prod_name, "coa") != NULL){
                for(int j = 0; j < sub_arr->size; j++){
                    char sub_name[NAME_SIZE];
                    strcpy(sub_name, metabolite_names[sub_arr->id[j]]);
                    remove_last_three_chars(sub_name);
                    if(strstr(sub_name, "accoa") != NULL){
                        remove_element(prod_arr, i);
                        remove_element(sub_arr, j);
                        return 0;
                    }
                }
            }
        }
        for(int i = 0; i < sub_arr->size; i++){
            char sub_name[NAME_SIZE];
            strcpy(sub_name, metabolite_names[sub_arr->id[i]]);
            remove_last_three_chars(sub_name);
            if(strstr(sub_name, "coa") !=NULL){
                for(int j = 0; j < prod_arr->size; j++){
                    char prod_name[NAME_SIZE];
                    strcpy(prod_name, metabolite_names[prod_arr->id[j]]);
                    remove_last_three_chars(prod_name);
                    if(strstr(prod_name, "malcoa") != NULL){
                        remove_element(sub_arr, i);
                        remove_element(prod_arr, j);
                        return 0;
                    }
                }
            }
        }
        for(int i = 0; i < prod_arr->size; i++){
            char prod_name[NAME_SIZE];
            strcpy(prod_name, metabolite_names[prod_arr->id[i]]);
            remove_last_three_chars(prod_name);
            if(strstr(prod_name, "coa") != NULL){
                for(int j = 0; j < sub_arr->size; j++){
                    char sub_name[NAME_SIZE];
                    strcpy(sub_name, metabolite_names[sub_arr->id[j]]);
                    remove_last_three_chars(sub_name);
                    if(strstr(sub_name, "malcoa") != NULL){
                        remove_element(prod_arr, i);
                        remove_element(sub_arr, j);
                        return 0;
                    }
                }
            }
        }
    }
    return 0;
}
void remove_last_three_chars(char *name) {
    int len = strlen(name); 
    if (len > 3) {
        name[len - 3] = '\0'; 
    } else {
        name[0] = '\0'; 
    }
}
int manual_delete_output(int reactnum,FILE*out,FILE* outnum,DynamicArray *sub_arr, DynamicArray *prod_arr, char **metabolite_names,char **reaction_names){
    printf("Manual deletion of cofactors:\n");
    while(sub_arr->size > 1 || prod_arr->size > 1){
        printf("For reaction: %s,we have more than one substrate and product.\n", reaction_names[reactnum]);
        printf("Substrates:\n");
        for(int i = 0; i < sub_arr->size; i++){
            char name[NAME_SIZE];
            strcpy(name, metabolite_names[sub_arr->id[i]]);
            printf("%d.%s\n",i,name);
        }
        printf("Products:\n");
        for(int i = 0; i < prod_arr->size; i++){
            char name[NAME_SIZE];
            strcpy(name, metabolite_names[prod_arr->id[i]]);
            printf("%d.%s\n",i+sub_arr->size,name);
        }
        printf("Please select the index of the metabolite to delete (0-%d),or input -1 to skip,or -2 to seperate: ", sub_arr->size + prod_arr->size - 1);
        int index;
        scanf("%d", &index);
        if (index == -1) {
            break; 
        } else if (index >= 0 && index < sub_arr->size) {
            remove_element(sub_arr, index);
        } else if (index >= sub_arr->size && index < sub_arr->size + prod_arr->size) {
            remove_element(prod_arr, index - sub_arr->size);
        } else if (index == -2){
            int sub_id, prod_id;
            double sub_value, prod_value;
            printf("For reaction: %s\n", reaction_names[reactnum]);
            printf("Substrates:\n");
            for(int i = 0; i < sub_arr->size; i++){
                char name[NAME_SIZE];
                strcpy(name, metabolite_names[sub_arr->id[i]]);
                printf("%d.%s\n",i,name);
            }
            printf("Products:\n");
            for(int i = 0; i < prod_arr->size; i++){
                char name[NAME_SIZE];
                strcpy(name, metabolite_names[prod_arr->id[i]]);
                printf("%d.%s\n",i,name);
            }
            printf("Please select the index of the substrate to seperate (0-%d): ", sub_arr->size - 1);
            scanf("%d", &sub_id);
            printf("Please select the index of the product to seperate (0-%d): ", prod_arr->size - 1);
            scanf("%d", &prod_id);
            for(int i = 0; i < sub_arr->size; i++){
                if(i == sub_id){
                    sub_id = sub_arr->id[i];
                    sub_value = sub_arr->rxn_data[i];
                    remove_element(sub_arr, i);
                    break;
                }
            }
            for(int i = 0; i < prod_arr->size; i++){
                if(i == prod_id){
                    prod_id = prod_arr->id[i];
                    prod_value = prod_arr->rxn_data[i];
                    remove_element(prod_arr, i);
                    break;
                }
            }
            DynamicArray *tsub_array = create_array(1);
            DynamicArray *tprod_array = create_array(1);
            add_element(tsub_array, sub_id, sub_value);
            add_element(tprod_array, prod_id, prod_value);
            output_reaction(reactnum, out, outnum, reaction_names, metabolite_names, tsub_array, tprod_array);
        }else{
            printf("Invalid index. Please try again.\n");
        }
    }
    return 0;
}
int output_reaction(int reactnum,FILE*out,FILE* outnum, char **reaction_names, char **metabolite_names,DynamicArray *sub_arr, DynamicArray *prod_arr){
    /*Reaction dataframe
    id m n
    a1 v1
    a2 v2
    ...
    am vm
    b1 w1
    b2 w2
    ...
    bn wn
    id: reaction id
    m: number of substrates
    n: number of products
    */
    fprintf(out, "%d\t%d\t%d\n", reactnum, sub_arr->size, prod_arr->size);
    for(int i = 0; i < sub_arr->size; i++){
        fprintf(out, "%d\t%lf\n", sub_arr->id[i], sub_arr->rxn_data[i]);
    }
    for(int i = 0; i < prod_arr->size; i++){
        fprintf(out, "%d\t%lf\n", prod_arr->id[i], prod_arr->rxn_data[i]);
    }
    fprintf(outnum, "%s\t%d\t%d\n", reaction_names[reactnum], sub_arr->size, prod_arr->size);
}
int ends_with_coa(const char *str) {
    size_t len = strlen(str);
    if (len >= 3 && strcmp(str + len - 3, "coa") == 0) {
        return 1; 
    }
    return 0; 
}