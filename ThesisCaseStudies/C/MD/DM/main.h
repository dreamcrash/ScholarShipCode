/* 
 * File:   main.h
 * Author: Dreamcrash
 *
 */

#ifndef MAIN_H
#define	MAIN_H

#ifdef	__cplusplus
extern "C" {
#endif

/** To validate the simulation */
void validate(MD *md){

     if(md->size < 3)
     {
        double refval[]={1731.4306625334357, 7397.392307839352, 16566.17177236039};
        double dev = fabs(md->ek-refval[md->size]);
        int interactions[] = {2506202, 53903910, 260577346};

        printf("[Ek,Dev] = [%f, %f]\n",md->ek,dev);
        printf("Inter => {%d}\n",md->interactions);

        if(dev >1.0e-10) printf("Simulation = {Fail!}\n");
        else
        {
            if(interactions[md->size] == md->interactions << 1)
            {
                printf("Simulation = {Valid!}\n");
            }
            else
            {
                printf("Fail\n");
            }
        }
     }
     else
     {
        printf("This size was not validated\n");
        printf("[Ek] = [%f]\n",md->ek);
     }
} 
   
int main                                (int argc, char** argv);
#ifdef	__cplusplus
}
#endif

#endif	/* MAIN_H */

