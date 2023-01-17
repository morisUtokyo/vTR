#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include "ksw2.h"

char *nextReadID;
int read_cnt;
#define BLK 4096
#define MAX_READ_LENGTH         100000
#define MAX_ID_LENGTH           10000
#define INTERNAL_UNIT_LENGTH    500
#define MAX_NUM_UNITS           10000

typedef struct{
    int     len;
    char    string[MAX_READ_LENGTH];
    char    info[MAX_READ_LENGTH];
    char    ID[MAX_READ_LENGTH];
    char    RegExpression[MAX_READ_LENGTH];
    char    pattern_string[MAX_READ_LENGTH];
    char    preciseRegExp[MAX_READ_LENGTH];
    int     numMatches;
    int     numMismatches;
    int     numDeletions;
    int     numInsertions;
    float   error_rate;
} Read;

typedef struct{
    int     len;
    char    type;
} cigarType;

typedef struct{
    char    unit[INTERNAL_UNIT_LENGTH];
    char    variant[INTERNAL_UNIT_LENGTH];
    int     occ;
    int     begin;
    int     posVariant;
    int     runLength;
    int     variantState;
} unitType;

// We obtained this function from https://github.com/lh3/ksw2/blob/master/ksw2.h and modified.
int align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape, cigarType *cigarArray)
{
    int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
    
    for (i = 0; i < ez.n_cigar; ++i){
        cigarArray[i].len   = ez.cigar[i]>>4;
        cigarArray[i].type  = "MID"[ez.cigar[i]&0xf];
        //printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
    }
    // Print the statistics of the alignment
    //printf("max_q=%d max_t=%d mqe=%d mqe_t=%d mte=%d mte_q=%d score=%d m_cigar=%d n_cigar=%d reach_end=%d\n", ez.max_q, ez.max_t, ez.mqe, ez.mqe_t, ez.mte, ez.mte_q, ez.score, ez.m_cigar, ez.n_cigar, ez.reach_end);
    free(ez.cigar); free(ts); free(qs);
    
    return(ez.n_cigar);
}


char capitalize(char c){
    char charCode;
    switch(c){
        case 'A':
        case 'a':
            charCode = 'A'; break;
        case 'C':
        case 'c':
            charCode = 'C'; break;
        case 'G':
        case 'g':
            charCode = 'G'; break;
        case 'T':
        case 't':
            charCode = 'T'; break;
        default:
            fprintf(stderr, "Invalid character: %c in capitalize\n", c); exit(EXIT_FAILURE);
    }
    return(charCode);
}

int pattern2string(char *a_pat, char *pattern_string, unitType *unitArray){
    char *pat;
    pat = (char *) malloc( sizeof(char) * strlen(a_pat) );
    strcpy(pat, a_pat);
    
    // pat=<AAG>12<AG>12 => AAG 12 AG 12
    for(int j=0; j<(int)strlen(pat); j++)
        if(pat[j]=='<' || pat[j]=='>') pat[j]=' ';
    
    char unit[100];
    int occ;
    int prevLen = strlen(pat);
    strcpy(pattern_string, ""); // Reset pattern_string
    
    int nUnits=0;
    int begin = 0;
    for(;;){
        // For example, retrieve "AAG 12" from " AAG 12 AG 12"
        sscanf(pat, " %s %d%[^\0]%*c", unit, &occ, pat);
        strcpy(unitArray[nUnits].unit, unit);
        unitArray[nUnits].occ = occ;
        unitArray[nUnits].begin = begin;
        unitArray[nUnits].posVariant = 0;
        unitArray[nUnits].runLength = 0;
        unitArray[nUnits].variantState = 0;
        nUnits++;
        begin += strlen(unit) * occ;
        
        // For example, concatenate 12 occurrences of AGG for AAG 12
        for(int i=0; i<occ; i++)
            sprintf(pattern_string, "%s%s", pattern_string, unit );
        int len = strlen(pat);
        if(prevLen == len) break; else prevLen = len;
    }
    free(pat);
    return(nUnits);
}

//#define DEBUG_unitVariant
int unitVariant(char cigarType, unitType *unitArray, Read *currentRead, int ip, int is, int iUnits, int *numMatches, int *numMismatches, int *numDeletions, int *numInsertions){
       
    char cp = currentRead->pattern_string[ip];
    char cs = currentRead->string[is];
    int  unitLen = (int) strlen(unitArray[iUnits].unit);
    int  begin = unitArray[iUnits].begin;
    int  occ = unitArray[iUnits].occ;
    int  posVariant = unitArray[iUnits].posVariant;
    
    if( posVariant >= INTERNAL_UNIT_LENGTH){
        fprintf(stderr, "The length of a unit variant %d exceeded %d.\n", posVariant, INTERNAL_UNIT_LENGTH);
        exit(1);
    }
    
    int new_ip = ip;
    switch(cigarType){
        case 'M':   // Increment ip, the index of the pattern.
                #ifdef DEBUG_unitVariant
                if( cp != cs  ) printf("* %d %c %c\n", (ip-begin)%unitLen, cp, cs);
                #endif
            unitArray[iUnits].variant[posVariant] = cs;
            if( cp != cs ){
                unitArray[iUnits].variantState=1;
                (*numMismatches)++;
            }else
                (*numMatches)++;
            posVariant++; new_ip++; break;
        case 'D':   // Increment ip.
                #ifdef DEBUG_unitVariant
                printf("* %d %c -\n", (ip-begin)%unitLen, cp);
                #endif
            unitArray[iUnits].variantState=1;
            (*numDeletions)++;
            new_ip++; break;
        case 'I':   // Do not increment ip.
                #ifdef DEBUG_unitVariant
                printf("* %d - %c\n", (ip-begin)%unitLen, cs);
                #endif
            unitArray[iUnits].variant[posVariant] = cs;
            unitArray[iUnits].variantState=1;
            (*numInsertions)++;
            posVariant++; break;
    }
    unitArray[iUnits].posVariant = posVariant;
    
    // Print the variant and reset posVariant to 0.
    if( 0 < ip && ((new_ip - begin)) % unitLen == 0 && cigarType != 'I' )
    {
        if(unitArray[iUnits].variantState == 0){
            (unitArray[iUnits].runLength)++;
        }else{
            // Print the previous run of representative units
            if(0 < unitArray[iUnits].runLength)
                sprintf( currentRead->preciseRegExp, "%s<%s>%d", currentRead->preciseRegExp, unitArray[iUnits].unit, unitArray[iUnits].runLength);
            unitArray[iUnits].runLength = 0;
            // Print the focal variant
            unitArray[iUnits].variant[posVariant] = '\0';
            sprintf( currentRead->preciseRegExp, "%s<%s>1", currentRead->preciseRegExp, unitArray[iUnits].variant);
        }
        // Reset
        unitArray[iUnits].variantState = 0;
        unitArray[iUnits].posVariant = 0;
    }
    
    if( begin + unitLen * occ <= new_ip ){
        // Print the previous run of representative units
        if(0 < unitArray[iUnits].runLength)
            sprintf( currentRead->preciseRegExp, "%s<%s>%d", currentRead->preciseRegExp, unitArray[iUnits].unit, unitArray[iUnits].runLength);
        // Move on to the next unit and reset run length etc.
        iUnits++;
        unitArray[iUnits].runLength = 0;
        unitArray[iUnits].variantState = 0;
        unitArray[iUnits].posVariant = 0;
    }
    return(iUnits);
}

//#define DEBUG_comp_preciseRegExp
//#define DEBUG_cigar
void comp_preciseRegExp(Read *currentRead){
    
    char *tmpS      = (char *) malloc(sizeof(char) * MAX_ID_LENGTH);
    char *content   = (char *) malloc(sizeof(char) * MAX_ID_LENGTH);
    char tag[100];
    
    strcpy(tmpS, currentRead->ID);
    #ifdef DEBUG_comp_preciseRegExp
    printf("comp_preciseRegExp\t%s\n", currentRead->ID);
    #endif
    int prev_tmpS_len= strlen(tmpS);
    int detected_a_pattern=0;
    for(;;){
        sscanf(tmpS, "#%s %s %[^\0]", tag, content, tmpS);
        #ifdef DEBUG_comp_preciseRegExp
        printf("tag=%s\tcontent=%s\ttmpS=%s\n", tag, content, tmpS);
        #endif
        if(strcmp(tag, "Info") == 0)
            strcpy(currentRead->info, content);
        if(strcmp(tag, "Pat") == 0){
            strcpy(currentRead->RegExpression, content);
            detected_a_pattern = 1;
        }
        if(prev_tmpS_len == (int)strlen(tmpS))  break;
        else    prev_tmpS_len = strlen(tmpS);
    }
    free(tmpS); free(content);
    if(detected_a_pattern == 0) return;
    
    for(int j=0; j<(int)strlen(currentRead->RegExpression); j++){
        if(currentRead->RegExpression[j]=='(') currentRead->RegExpression[j]='<';
        if(currentRead->RegExpression[j]==')') currentRead->RegExpression[j]='>';
    }
    
    strcpy(currentRead->preciseRegExp, "");
    cigarType *cigarArray = (cigarType *) malloc(sizeof(cigarType) * MAX_READ_LENGTH);
    unitType  *unitArray  = (unitType *) malloc(sizeof(unitArray) * MAX_NUM_UNITS);
    
    // Call ksw2 to calculate the alignment score
    int nUnits = pattern2string(currentRead->RegExpression, currentRead->pattern_string, unitArray);
    int len_pattern = 0;
    for(int j=0; j<nUnits; j++)
        len_pattern += strlen(unitArray[j].unit) * unitArray[j].occ;
    int sc_mch = 1; int sc_mis = -1; int gapo = 0; int gape = 1;
    int nCigar = align(currentRead->pattern_string, currentRead->string, sc_mch, sc_mis, gapo, gape, cigarArray);
        
    // While scanning the CIGAR format, calculate variants of each representative unit.
    int ip = 0;  // ip and is scan the pattern and given string.
    int is = 0;
    int iUnits=0;
    int numMatches, numMismatches, numDeletions, numInsertions;
    numMatches = numMismatches = numInsertions = numDeletions = 0;
    for(int j=0; j<nCigar; j++){
        if(len_pattern < ip+cigarArray[j].len || currentRead->len < is+cigarArray[j].len) break;
        char cigarType = cigarArray[j].type;
            #ifdef DEBUG_cigar
            printf("%c %d\n", cigarType, cigarArray[j].len);
            #endif
        switch(cigarType){
            case 'M':   // Increment ip, the index of the pattern.
                for(int k=0; k<cigarArray[j].len; k++)
                    iUnits = unitVariant(cigarType, unitArray, currentRead, ip+k, is+k, iUnits, &numMatches, &numMismatches, &numDeletions, &numInsertions);
                ip += cigarArray[j].len;
                is += cigarArray[j].len;
                break;
            case 'D':   // Increment ip.
                for(int k=0; k<cigarArray[j].len; k++)
                    iUnits = unitVariant(cigarType, unitArray, currentRead, ip+k, is, iUnits, &numMatches, &numMismatches, &numDeletions, &numInsertions);
                ip += cigarArray[j].len;
                break;
            case 'I':   // Do not increment ip.
                for(int k=0; k<cigarArray[j].len; k++)
                    iUnits = unitVariant(cigarType, unitArray, currentRead, ip, is+k, iUnits, &numMatches, &numMismatches, &numDeletions, &numInsertions);
                is += cigarArray[j].len;
                break;
        }
    }
    currentRead->numMatches     = numMatches;
    currentRead->numMismatches  = numMismatches;
    currentRead->numDeletions   = numDeletions;
    currentRead->numInsertions  = numInsertions;
    currentRead->error_rate     = (float) (numMismatches+numDeletions+numInsertions)/(numMatches+numMismatches+numDeletions+numInsertions);
    fflush(stdout);
    free(cigarArray);
    free(unitArray);
}


void return_one_read(FILE *fp, Read *currentRead){
    
    char *s = (char *) malloc( sizeof(char) * MAX_READ_LENGTH );
    int i;
    int cnt=0;
    int no_read = 1;
    
    while (fgets(s, BLK, fp) != NULL) { // Feed a string of size BLK from fp into string s
        fflush(fp);
        no_read = 0;
        
        if(s[0] == '>'){
            if(read_cnt != -1){ // This is NOT the first read
                // Set the ID of currentRead to the ID of nextRead
                int j;
                for(j=0;  nextReadID[j] != '\0'; j++)
                    currentRead->ID[j] = nextReadID[j];
                currentRead->ID[j] = '\0';
            }
            // Feed the ID of the current read into the ID of nextRead
            int shift;
            if(s[1]==' ') shift=2; else shift=1;  // Skip the space at the head
            for(i=0; s[shift+i]!='\0' && s[shift+i]!='\n' && s[shift+i]!='\r' && i<BLK; i++)
                nextReadID[i] = s[shift+i];
            nextReadID[i] = '\0';
            
            if(read_cnt == -1){ // This is the first read
                read_cnt = 0;
            }else{
                read_cnt++;
                // Finalize the currentRead string by appending '\0'
                currentRead->string[cnt] = '\0';
                currentRead->len = cnt;
                return;
            }
        }else{
            // Feed the string
            for(i=0; s[i]!='\0' && s[i]!='\n' && s[i]!='\r' && i<BLK; i++){
                currentRead->string[cnt] = capitalize(s[i]);
                cnt++;
                if( MAX_READ_LENGTH <= cnt){
                    fprintf(stderr, "fatal error: The length %d is tentatively at most %i.\nread ID = %s\nSet MAX_READ_LENGTH to a larger value", cnt, MAX_READ_LENGTH, currentRead->ID);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    if(no_read == 1){  // No reads
        currentRead->len = 0;
    }else{
        // Process the last read.
        // Set the ID of currentRead to the ID of nextRead
        int j;
        for(j=0;  nextReadID[j] != '\0'; j++)
            currentRead->ID[j] = nextReadID[j];
        currentRead->ID[j] = '\0';
        // Finalize the currentRead string by appending '\0'
        currentRead->string[cnt] = '\0';
        currentRead->len = cnt;
        read_cnt++;
    }
    free(s);
}



int main(int argc, char *argv[])
{
    char fastaFileName[1000];
    int opt;
    while ((opt = getopt(argc, argv, "f:")) != -1) {
        switch(opt){
            case 'f':   strcpy(fastaFileName,optarg); break;
            default:    exit(EXIT_FAILURE);
        }
    }
    FILE *fp = fopen(fastaFileName, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", fastaFileName);
        exit(EXIT_FAILURE);
    }

    // Feed one read and process it iteratively.
    Read *currentRead = (Read *) malloc(sizeof(Read));
    if(currentRead==NULL){
        fprintf(stderr, "Failure to malloc currentRead\n");
        exit(EXIT_FAILURE);
    }
    nextReadID = (char *) malloc( sizeof(char) * (MAX_ID_LENGTH+1) );
    if(nextReadID == NULL) exit(EXIT_FAILURE);
    read_cnt = -1;
    
    for(int i=0;;i++){
        return_one_read(fp, currentRead);
        if(currentRead->len == 0) break;

        comp_preciseRegExp(currentRead);
        printf("> #Len %d #Err %3.3f #Pat %s #PrecisePat %s\n%s\n", currentRead->len, currentRead->error_rate, currentRead->RegExpression, currentRead->preciseRegExp, currentRead->string);
    }
    free(currentRead);
    free(nextReadID);

    exit(0);
}
