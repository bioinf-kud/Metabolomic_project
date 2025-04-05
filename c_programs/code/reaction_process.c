#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
struct enzyme_info{//store enzyme data
    int enzyme_mass;
    char enzyme_gene_id_prot[100];
    char enzyme_gene_id_trans[100];
    char Uniprot_id[20];
    double km;
    double kcat;
    double enzyme_MS;
    double enzyme_conc;  
    double enzyme_effeciency;
};
struct enzyme{//scan enzyme data from file and store temporarily
    char reaction_id[100];
    char substrate_id[100];
    char enzyme_gene_id_prot[100];
    char enzyme_gene_id_trans[100];
    char Uniprot_id[20];
    int enzyme_mass;
    double km;
    double kcat; 
};
struct reaction_info{//store reaction data
    char reaction_id[100];
    char substrate_id[100];
    double substrate_conc;
    double enzyme_calc_value;
    int enzyme_num;
    struct enzyme_info *enzyme;
};
//Scan number of reaction in the file, reactions must be list in order, and provide related substrate information.
int scan_reaction_info(char *path, struct reaction_info *reaction_list);  
//Read reaction data from file and store in reaction_list.
int read_reaction_info(char *path, struct reaction_info *reaction_list,int scnt);
//Read a line from enzyme file.
struct enzyme read_enzyme_info(FILE *f);
//Scan the enzyme file to find the number of enzymes related to a specific reaction and substrate.
int scan_enzyme_list_info(char *path, char* reaction, char*substrate);
//Read enzyme data from file and store in reaction_list.
int read_enzyme_list(char *path, struct reaction_info *reaction_list);
//Validate the theory by proteomics data.
int calc_all_enz_MS(char *MS_path,char *out_path,struct reaction_info *reaction_list,int reaction_num);
//Read proteomics data from file and store in reaction_list.
int read_enz_MS(FILE *f,struct reaction_info *reaction_list,char**enzyme_gene_id_list,int enzcnt,int reaction_num);
//Calculate enzyme concentration and efficiency.
int calc_enz_MS(struct reaction_info *reaction_list,int reaction_num);
//Output the result to a file.
int output_enz_MS(FILE *f,struct reaction_info *reaction_list,char*name,int reaction_num);
//Validate the theory by transcriptom data.
int calc_all_enz_TS(char *TS_path,char *out_path,struct reaction_info *reaction_list,int reaction_num);
//Read proteomics data from file and store in reaction_list.
int read_enz_TS(FILE *f,struct reaction_info *reaction_list,char**enzyme_gene_id_list,int enzcnt,int reaction_num);
//Calculate enzyme concentration and efficiency.
int calc_enz_TS(struct reaction_info *reaction_list,int reaction_num);
int main(int argc, char *argv[]){
    char reaction_file_path[200];
    char enzyme_file_path[200];
    char enzyme_MS_file_path[200];
    char enzyme_MS_out_file_path[200];
    char enzyme_TS_file_path[200];
    char enzyme_TS_out_file_path[200];
    int flag=0,cnt=0,flag2=0;
    int opt;
    while ((opt = getopt(argc, argv, "r:e:i:o:hpt")) != -1) {//Get options from command line.
        switch (opt) {
            case 'r': 
                strcpy(reaction_file_path,optarg);
                cnt++;
                break;
            case 'e': 
                strcpy(enzyme_file_path,optarg);
                cnt++;
                break;
            case 'i': 
                strcpy(enzyme_MS_file_path,optarg);
                strcpy(enzyme_TS_file_path,optarg);
                cnt++;
                break;
            case 'o': 
                strcpy(enzyme_MS_out_file_path,optarg);
                strcpy(enzyme_TS_out_file_path,optarg);
                cnt++;
                break;
            case 'p': 
                flag=1;
                flag2+=1;
                break;
            case 't':
                flag=2;
                flag2+=1;
                break;
            case 'h': 
                printf("Usage: %s [-r reaction_file_path] [-e enzyme_information_file_path] [-i enzyme_omics_data_file_path] [-o output_file_path] [-p] [-t]\n", argv[0]);
                printf("Options:\n");
                printf("  -r reaction_file_path: Path to the reaction file.\n");
                printf("  -e enzyme_information_file_path: Path to the enzyme file.\n");
                printf("  -i enzyme_omics_data_file_path: Path to the proteomics or transcriptomics data file.\n");
                printf("  -o output_file_path: Path to the output file.\n");
                printf("  -p: Use proteomics data to validate the theory.\n");
                printf("  -t: Use transcriptom data to validate the theory.\n");
                printf("  -h: Display this help message.\n");
                return 0;
            case '?': // 未知选项
                printf("Unknown option: %c\n", optopt);
                break;
            default:
            printf("Usage: %s [-r reaction_file_path] [-e enzyme_information_file_path] [-i enzyme_omics_data_file_path] [-o output_file_path] [-p] [-t]\n", argv[0]);
                return 1;
        }
    }
    if(cnt!=4){//Check if all necessary arguments are provided.
        printf("Error:Not enough arguments!\n");
        return 0;
    }
    if(flag2!=1){//Check if only one omics method is selected.
        printf("Error:No omics method or too many selected!\n");
        return 0;
    }
    struct reaction_info *reaction_list;//List of reactions in a reaction chain of interest.
    printf("Reading reaction data...\n");
    int reaction_num=scan_reaction_info(reaction_file_path,reaction_list);
    reaction_list = (struct reaction_info*)malloc(sizeof(struct reaction_info)*reaction_num);
    read_reaction_info(reaction_file_path,reaction_list,reaction_num);
    printf("Reading enzyme data...\n");
    for(int i=0;i<reaction_num;i++){
        int enzyme_num = scan_enzyme_list_info(enzyme_file_path,reaction_list[i].reaction_id,reaction_list[i].substrate_id);
        reaction_list[i].enzyme = (struct enzyme_info*)malloc(sizeof(struct enzyme_info)*enzyme_num);
        reaction_list[i].enzyme_num = enzyme_num;
        if(enzyme_num==0){
            printf("Error:No enzyme data found for reaction %s.\n",reaction_list[i].reaction_id);
            return 0;
        }
    }
    for(int i=0;i<reaction_num;i++){
        int a=read_enzyme_list(enzyme_file_path,reaction_list+i);
        if(a==1)
            return 0;
    }
    printf("Calculating data...\n");
    int ctr;
    switch(flag){//Select the omics method.
        case 1:
            printf("Using proteomics data...\n");
            ctr=calc_all_enz_MS(enzyme_MS_file_path,enzyme_MS_out_file_path,reaction_list,reaction_num);
            break;
        case 2:
            printf("Using transcriptom data...\n");
            ctr=calc_all_enz_TS(enzyme_TS_file_path,enzyme_TS_out_file_path,reaction_list,reaction_num);
            break;
        default:
            printf("Error:No omics data selected!\n");
            return 0;
    }
    if (ctr==1)
        return 0;
    printf("Data calculated successfully!\n");
    return 0;
}
int scan_reaction_info(char *path, struct reaction_info *reaction_list){
    FILE *scan;
    scan = fopen(path,"r");
    if(scan==NULL){
        printf("Reaction file not found\n");
        return 1;
    }
    int scnt=0;
    while(1){
        char c=fgetc(scan);
        if(c==EOF)
            break;
        if(c=='\n')
            scnt++;
    }
    fclose(scan);
    return scnt;
}
int read_reaction_info(char *path, struct reaction_info *reaction_list,int scnt){
    FILE *read;
    read = fopen(path,"r");
    char temp_reaction_id[100];
    char temp_substrate_id[100];
    double temp_substrate_conc;
    double temp_Kadj;
    for(int i=0;i<scnt;i++){
        int recnt=0;
        while(1){
            char c=fgetc(read);
            if(c==EOF){
                if(i!=scnt-1){
                    printf("Error reading reaction!\n");
                    return 0;
                }
                return scnt;
            }
            temp_reaction_id[recnt]=c;
            if(c=='\t'){
                temp_reaction_id[recnt]='\0';
                break;
            }
            recnt++;
        }
        int sucnt=0;
        while(1){
            char c=fgetc(read);
            temp_substrate_id[sucnt]=c;
            if(c=='\t'){
                temp_substrate_id[sucnt]='\0';
                break;
            }
            sucnt++;
        }
        fscanf(read,"%lf\n",&temp_substrate_conc);
        strcpy(reaction_list[i].reaction_id,temp_reaction_id);
        strcpy(reaction_list[i].substrate_id,temp_substrate_id);
        reaction_list[i].substrate_conc = temp_substrate_conc;
    }
    return scnt;
}
struct enzyme read_enzyme_info(FILE *f){
    struct enzyme temp_enzyme;
    int recnt=0;
    while(1){
        char c=fgetc(f);
        if(c==EOF){
            strcpy(temp_enzyme.reaction_id,"EOFEOF");
            return temp_enzyme;
        }
        temp_enzyme.reaction_id[recnt]=c;
        if(c=='\t'){
            temp_enzyme.reaction_id[recnt]='\0';
            break;
        }
        recnt++;
    }
    int sucnt=0;
    while(1){
        char c=fgetc(f);
        temp_enzyme.substrate_id[sucnt]=c;
        if(c=='\t'){
            temp_enzyme.substrate_id[sucnt]='\0';
            break;
        }
        sucnt++;
    }
    int gepcnt=0;
    while(1){
        char c=fgetc(f);
        temp_enzyme.enzyme_gene_id_prot[gepcnt]=c;
        if(c=='\t'){
            temp_enzyme.enzyme_gene_id_prot[gepcnt]='\0';
            break;
        }
        gepcnt++;
    }
    int gecnt=0;
    while(1){
        char c=fgetc(f);
        temp_enzyme.enzyme_gene_id_trans[gecnt]=c;
        if(c=='\t'){
            temp_enzyme.enzyme_gene_id_trans[gecnt]='\0';
            break;
        }
        gecnt++;
    }
    int upcnt=0;
    while(1){
        char c=fgetc(f);
        temp_enzyme.Uniprot_id[upcnt]=c;
        if(c=='\t'){
            temp_enzyme.Uniprot_id[upcnt]='\0';
            break;
        }
        upcnt++;
    }
        fscanf(f,"%d\t%lf\t%lf\n",&temp_enzyme.enzyme_mass,&temp_enzyme.km,&temp_enzyme.kcat);
    return temp_enzyme;
}
int scan_enzyme_list_info(char *path, char* reaction, char*substrate){
    FILE *scan;
    scan = fopen(path,"r");
    if(scan==NULL){
        printf("Enzyme file not found!\n");
        return 0;
    }
    int scnt=0;
    struct enzyme temp_enzyme;
    while(1){
        temp_enzyme = read_enzyme_info(scan);
        if(strcmp(temp_enzyme.reaction_id,"EOFEOF")==0)
            break;
        if(strcmp(temp_enzyme.reaction_id,reaction)==0 && strcmp(temp_enzyme.substrate_id,substrate)==0)
            scnt++;
    }
    fclose(scan);
    return scnt;
}
int read_enzyme_list(char *path, struct reaction_info *reaction_list){
    FILE *read;
    read = fopen(path,"r");
    struct enzyme temp_enzyme;
    for(int i=0;i<reaction_list->enzyme_num;i++){
        while(1){
            temp_enzyme = read_enzyme_info(read);
            if(strcmp(temp_enzyme.reaction_id,"EOFEOF")==0){
                if(i!=reaction_list->enzyme_num-1){
                    printf("Error reading enzyme!\n");
                    return 1;
                }
                return 0;
            }
            if(strcmp(temp_enzyme.reaction_id,reaction_list->reaction_id)==0 && strcmp(temp_enzyme.substrate_id,reaction_list->substrate_id)==0){
                strcpy(reaction_list->enzyme[i].enzyme_gene_id_prot,temp_enzyme.enzyme_gene_id_prot);
                strcpy(reaction_list->enzyme[i].enzyme_gene_id_trans,temp_enzyme.enzyme_gene_id_trans);
                strcpy(reaction_list->enzyme[i].Uniprot_id,temp_enzyme.Uniprot_id);
                reaction_list->enzyme[i].enzyme_mass = temp_enzyme.enzyme_mass;
                reaction_list->enzyme[i].km = temp_enzyme.km;
                reaction_list->enzyme[i].kcat = temp_enzyme.kcat;
                break;
            }    
        } 
    }
    fclose(read);
    return 0;
}
int calc_all_enz_MS(char *MS_path,char *out_path,struct reaction_info *reaction_list,int reaction_num){
    //Scan proteomics data file to find the number of enzymes and samples.
    FILE *scanMS,*readMS;
    scanMS = fopen(MS_path,"r");
    if(scanMS==NULL){
        printf("Proteomics data not found!\n");
        return 1;
    }
    int enzcnt=0;
    int sampcnt=1;//May cause error if the last line is empty. To debug, change to 0.
    char **enzyme_gene_id_list;
    while(1){
        char c=fgetc(scanMS);
        if(c=='\t'){
            enzcnt++;
        }        
        if(c=='\n'){
            break;
        }
    }
    while(1){
        char c=fgetc(scanMS);
        if(c=='\n'){
            sampcnt++;
        }        
        if(c==EOF){
            break;
        }
    }
    printf("%d\n",enzcnt);
    fclose(scanMS);
    readMS = fopen(MS_path,"r");
    //read and store all enzyme names.
    enzyme_gene_id_list = (char**)malloc(sizeof(char*)*enzcnt);
    for(int i=0;i<enzcnt;i++){
        enzyme_gene_id_list[i] = (char*)malloc(sizeof(char)*100);
    }
    while(1){
        char c=fgetc(readMS);
        if(c=='\t'){
            break;
        }
    }
    for(int i=0;i<enzcnt;i++){
        int cnt=0;
        while(1){
            char c=fgetc(readMS);
            if(c=='\t' || c=='\n'){
                enzyme_gene_id_list[i][cnt]='\0';
                break;
            }
            enzyme_gene_id_list[i][cnt]=c;
            cnt++;
        }
        printf("%d",i);
    }
    FILE *out;
    out = fopen(out_path,"w");
    //Read line by line and calculate the enzyme concentration and efficiency of each sample.
    for(int i=0;i<sampcnt;i++){
        char name[100];
        int cnt=0;
        while(1){
            char c=fgetc(readMS);
            if(c=='\t'){
                name[cnt]='\0';
                break;
            }
            name[cnt]=c;
            cnt++;
        }
        printf("Processing:%s.\n",name);
        int a=read_enz_MS(readMS,reaction_list,enzyme_gene_id_list,enzcnt,reaction_num);
        fgetc(readMS);
        fgetc(readMS);
        if(a==1)
            return 1;
        calc_enz_MS(reaction_list,reaction_num);  
        //Output the result to a file.
        if(strstr(name,"bridge") == NULL)
            output_enz_MS(out,reaction_list,name,reaction_num);
        printf("\n");  
    }
    fclose(readMS);
    fclose(out);
    return 0;
}
int read_enz_MS(FILE *f,struct reaction_info *reaction_list,char**enzyme_gene_id_list,int enzcnt,int reaction_num){
    for(int j=0;j<reaction_num;j++)
        for(int k=0;k<reaction_list[j].enzyme_num;k++)
            reaction_list[j].enzyme[k].enzyme_MS =0;
    for(int i=0;i<enzcnt;i++){
        double temp_MS;
        fscanf(f,"%lf",&temp_MS);
        int rcnt=0,ecnt=0;
        for(int j=0;j<reaction_num;j++)
            for(int k=0;k<reaction_list[j].enzyme_num;k++)
                if(strcmp(reaction_list[j].enzyme[k].enzyme_gene_id_prot,enzyme_gene_id_list[i])==0){
                    reaction_list[j].enzyme[k].enzyme_MS += temp_MS;
                    break;
                }  
    }
    int rcnt=0;
    for(int i=0;i<reaction_num;i++){
        rcnt=0;
        for(int j=0;j<reaction_list[i].enzyme_num;j++){
            if(reaction_list[i].enzyme[j].enzyme_MS!=0){
                rcnt++;
            }
            else{
                printf("Warning:No proteomics data found for enzyme %s in reaction %s.\n",reaction_list[i].enzyme[j].Uniprot_id,reaction_list[i].reaction_id);
            }
        }
        if(rcnt==0){
            printf("Error:No proteomics data found for reaction %s.\n",reaction_list[i].reaction_id);
            return 1;
        }
    }
    return 0;
}
int calc_enz_MS(struct reaction_info *reaction_list,int reaction_num){
    for(int i=0;i<reaction_num;i++){
        for(int j=0;j<reaction_list[i].enzyme_num;j++){
            reaction_list[i].enzyme[j].enzyme_conc = reaction_list[i].enzyme[j].enzyme_MS/reaction_list[i].enzyme[j].enzyme_mass;
            reaction_list[i].enzyme[j].enzyme_effeciency = reaction_list[i].enzyme[j].kcat*1000/(reaction_list[i].enzyme[j].km+reaction_list[i].substrate_conc);
        }
    }
    for(int i=0;i<reaction_num;i++){
        double temp_e=0,temp_k=0;
        for(int j=0;j<reaction_list[i].enzyme_num;j++){
            temp_e += reaction_list[i].enzyme[j].enzyme_conc;
            temp_k += reaction_list[i].enzyme[j].enzyme_effeciency*reaction_list[i].enzyme[j].enzyme_conc;
        }
        reaction_list[i].enzyme_calc_value=temp_e*temp_k;
    }
    return 0;
}
int output_enz_MS(FILE *f,struct reaction_info *reaction_list,char*name,int reaction_num){
    fprintf(f,"%s\t",name);
    for(int i=0;i<reaction_num-1;i++){
        fprintf(f,"%.10lf\t",reaction_list[i].enzyme_calc_value/reaction_list[i+1].enzyme_calc_value);
    }
    fprintf(f,"\n");
    return 0;
}
int calc_all_enz_TS(char *TS_path,char *out_path,struct reaction_info *reaction_list,int reaction_num){
    FILE *scanTS,*readTS;
    scanTS = fopen(TS_path,"r");
    if(scanTS==NULL){
        printf("Transcriptom data not found!\n");
        return 1;
    }
    int enzcnt=0;
    int sampcnt=0;//May cause error if the last line is empty.
    char **enzyme_gene_id_list;
    while(1){
        char c=fgetc(scanTS);
        if(c=='\t'){
            enzcnt++;
        }        
        if(c=='\n'){
            break;
        }
    }
    while(1){
        char c=fgetc(scanTS);
        if(c=='\n'){
            sampcnt++;
        }        
        if(c==EOF){
            break;
        }
    }
    printf("%d\n",enzcnt);
    fclose(scanTS);
    readTS = fopen(TS_path,"r");
    enzyme_gene_id_list = (char**)malloc(sizeof(char*)*enzcnt);
    for(int i=0;i<enzcnt;i++){
        enzyme_gene_id_list[i] = (char*)malloc(sizeof(char)*100);
    }
    while(1){
        char c=fgetc(readTS);
        if(c=='\t'){
            break;
        }
    }
    for(int i=0;i<enzcnt;i++){
        int cnt=0;
        while(1){
            char c=fgetc(readTS);
            if(c=='\t' || c=='\n'){
                enzyme_gene_id_list[i][cnt]='\0';
                break;
            }
            enzyme_gene_id_list[i][cnt]=c;
            cnt++;
        }
    }
    for(int i=0;i<enzcnt;i++){
        printf("%s\n",enzyme_gene_id_list[i]);
    }
    FILE *out;
    out = fopen(out_path,"w");
    for(int i=0;i<sampcnt;i++){
        char name[100];
        int cnt=0;
        while(1){
            char c=fgetc(readTS);
            if(c=='\t'){
                name[cnt]='\0';
                break;
            }
            name[cnt]=c;
            cnt++;
        }
        printf("Processing:%s.\n",name);
        int a=read_enz_TS(readTS,reaction_list,enzyme_gene_id_list,enzcnt,reaction_num);
        fgetc(readTS);
        if(a==1)
            return 1;
        calc_enz_TS(reaction_list,reaction_num);  
        if(strstr(name,"bridge") == NULL)
            output_enz_MS(out,reaction_list,name,reaction_num);
        printf("\n");   
    }
    fclose(readTS);
    fclose(out);
    return 0;
}
int read_enz_TS(FILE *f,struct reaction_info *reaction_list,char**enzyme_gene_id_list,int enzcnt,int reaction_num){
    for(int j=0;j<reaction_num;j++)
        for(int k=0;k<reaction_list[j].enzyme_num;k++)
            reaction_list[j].enzyme[k].enzyme_MS =0;
    for(int i=0;i<enzcnt;i++){
        double temp_MS;
        fscanf(f,"%lf",&temp_MS);
        int rcnt=0,ecnt=0;
        for(int j=0;j<reaction_num;j++)
            for(int k=0;k<reaction_list[j].enzyme_num;k++)
                if(strcmp(reaction_list[j].enzyme[k].enzyme_gene_id_trans,enzyme_gene_id_list[i])==0){
                    reaction_list[j].enzyme[k].enzyme_MS += temp_MS;
                    break;
                }  
    }
    int rcnt=0;
    for(int i=0;i<reaction_num;i++){
        rcnt=0;
        for(int j=0;j<reaction_list[i].enzyme_num;j++){
            if(reaction_list[i].enzyme[j].enzyme_MS!=0){
                rcnt++;
            }
            else{
                printf("Warning:No transcriptom data found for enzyme %s in reaction %s.\n",reaction_list[i].enzyme[j].Uniprot_id,reaction_list[i].reaction_id);
            }
        }
        if(rcnt==0){
            printf("Error:No transcriptom data found for reaction %s.\n",reaction_list[i].reaction_id);
            return 1;
        }
    }
    return 0;
}
int calc_enz_TS(struct reaction_info *reaction_list,int reaction_num){
    for(int i=0;i<reaction_num;i++){
        for(int j=0;j<reaction_list[i].enzyme_num;j++){
            reaction_list[i].enzyme[j].enzyme_conc = reaction_list[i].enzyme[j].enzyme_MS;
            reaction_list[i].enzyme[j].enzyme_effeciency = reaction_list[i].enzyme[j].kcat*1000/(reaction_list[i].enzyme[j].km+reaction_list[i].substrate_conc);
        }
    }
    for(int i=0;i<reaction_num;i++){
        double temp_e=0,temp_k=0;
        for(int j=0;j<reaction_list[i].enzyme_num;j++){
            temp_e += reaction_list[i].enzyme[j].enzyme_conc;
            temp_k += reaction_list[i].enzyme[j].enzyme_effeciency*reaction_list[i].enzyme[j].enzyme_conc;
        }
        reaction_list[i].enzyme_calc_value=temp_e*temp_k;
    }
    return 0;
}
