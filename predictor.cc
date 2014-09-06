#include "predictor.h"
#include <cstdlib>
#include <time.h>
#include <bitset>

#define BIMODAL_CTR_MAX  3
#define BIMODAL_CTR_INIT 2
#define TAGPRED_CTR_MAX  7
#define TAGPRED_CTR_INIT 0
#define BIMODALLOG   14
#define NUMTAGTABLES 4
#define TAGPREDLOG 12

/////////////// STORAGE BUDGET JUSTIFICATION ////////////////
// Total storage budget: 32KB + 17 bits
// Total PHT counters: 2^17 
// Total PHT size = 2^17 * 2 bits/counter = 2^18 bits = 32KB
// GHR size: 17 bits
// Total Size = PHT size + GHR size
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PREDICTOR::PREDICTOR(void)
{
    
 /// First initiating Bimodal Table
   // Its a simple 2 bit counter table 
    
  bimodalLog    = BIMODALLOG;
  numBimodalEntries   = (1<< bimodalLog);
    bimodal = new UINT32[numBimodalEntries];

     for(UINT32 ii=0; ii< numBimodalEntries; ii++)
     {
        bimodal[ii]=BIMODAL_CTR_INIT;
     }
    
   // Next to initiating the taggedPredictors
    tagPredLog = TAGPREDLOG;
    numTagPredEntries = (1 << tagPredLog);
    //cout << " No of entries in tag predictors = " << numTagPredEntries << endl;          
    
    for(UINT32 ii = 0; ii < NUMTAGTABLES ; ii++)
    {
       tagPred[ii] = new TagEntry[numTagPredEntries];
    }
    for(UINT32 ii = 0; ii < NUMTAGTABLES; ii++)
    {
        for(UINT32 j =0; j < numTagPredEntries; j++)
        {
            tagPred[ii][j].ctr = 0;
            tagPred[ii][j].tag = 0;
            tagPred[ii][j].usefulBits = 0;
        }
      
    }
    
    // Geometric lengths of history taken to consider correlation of different age.
    // Table 0 with the longest history as per PPM code
    geometric[0] = 130;
    geometric[1] = 44;
    geometric[2] = 15;
    geometric[3] = 5;
    /*this gives 3.41 MPKI !!
     * geometric[0] = 200;
    geometric[1] = 80;
    geometric[2] = 20;
    geometric[3] = 5;*/ 
    
    // Initializing Compressed Buffers.
    // first for index of the the tagged tables
    for(int i = 0; i < NUMTAGTABLES; i++)
    {
        indexComp[i].compHist = 0;
        indexComp[i].geomLength = geometric[i]; 
        indexComp[i].targetLength = TAGPREDLOG;
    }
    
    // next for the tag themselves
    
        // The tables have different tag lengths
        // T2 and T3 have tag length -> 8
        // T0 and T1 have tag length -> 9
        // second index indicates the Bank no.
        // Reason for using two -> PPM paper... single 
  // CSR is sensitive to periodic patterns in global history which is a common case
    for(int j = 0; j < 2 ; j++)
    {
        for(int i=0; i < NUMTAGTABLES; i++)
        {
            tagComp[j][i].compHist = 0;
            tagComp[j][i].geomLength = geometric[i];
            if(j == 0)
            {
                if(i < 2)
                tagComp[j][i].targetLength = 9 ;
                else
                tagComp[j][i].targetLength = 9 ;    
            }
            else
            {
                if(i < 2)
                tagComp[j][i].targetLength = 8 ;
                else
                tagComp[j][i].targetLength = 8 ;
            }
        }   
    }    
       // Preditions banks and prediction values 
       primePred = -1;
       altPred = -1;
       primeBank = NUMTAGTABLES;
       altBank = NUMTAGTABLES;
       
       for(int i=0; i < NUMTAGTABLES; i++)
       {    
            indexTagPred[i] = 0;
       }
       for(int i=0; i < NUMTAGTABLES; i++)
       {    
            tag[i] = 0;
       }
       clock = 0;
       clock_flip = 1;
       PHR = 0;
       GHR.reset();
       altBetterCount = 8;
}      

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

bool   PREDICTOR::GetPrediction(UINT32 PC){
  
 /// Base Prediction   
    
  bool basePrediction;
  UINT32 bimodalIndex   = (PC) % (numBimodalEntries);
  UINT32 bimodalCounter = bimodal[bimodalIndex];
  
  if(bimodalCounter > BIMODAL_CTR_MAX/2){
    basePrediction =  1; 
  }else{
    basePrediction =  0; 
  }
  
      // Hash to get tag includes info about bank, pc and global history compressed
    
    // formula given in PPM paper 
    // pc[9:0] xor CSR1 xor (CSR2 << 1)
    for(int i = 0; i < NUMTAGTABLES; i++)
    { 
       tag[i] = PC ^ tagComp[0][i].compHist ^ (tagComp[1][i].compHist << 1);
       /// These need to be masked
     // 9 bit tags for T0 and T1 // 8 bit tags for T2 and T3
    }
    tag[0] &= ((1<<9)-1);
    tag[1] &= ((1<<9)-1);
    tag[2] &= ((1<<9)-1);
    tag[3] &= ((1<<9)-1);
    // Tags now ready
    
  // How to get index for each bank ??
  // bank 1
    // Hash of PC, PC >> index length important , GHR geometric, path info
    UINT32 index_mask = ((1<<TAGPREDLOG) - 1);
    
           
       /*indexTagPred[0] = PC ^ (PC >> TAGPREDLOG) ^ indexComp[0].compHist ^ PHR ^ (PHR >> TAGPREDLOG);
       indexTagPred[1] = PC ^ (PC >> TAGPREDLOG) ^ indexComp[1].compHist ^ PHR ^ (PHR >> TAGPREDLOG);
       indexTagPred[2] = PC ^ (PC >> TAGPREDLOG) ^ indexComp[2].compHist ^ PHR ^ (PHR >> TAGPREDLOG);
       indexTagPred[3] = PC ^ (PC >> TAGPREDLOG) ^ indexComp[3].compHist ^ PHR ^ (PHR >> TAGPREDLOG);*/
       indexTagPred[0] = PC ^ (PC >> TAGPREDLOG) ^ indexComp[0].compHist ^ PHR ^ (PHR >> TAGPREDLOG);
       indexTagPred[1] = PC ^ (PC >> (TAGPREDLOG - 1)) ^ indexComp[1].compHist ^ (PHR );
       indexTagPred[2] = PC ^ (PC >> (TAGPREDLOG - 2)) ^ indexComp[2].compHist ^ (PHR & 31);
       indexTagPred[3] = PC ^ (PC >> (TAGPREDLOG - 3)) ^ indexComp[3].compHist ^ (PHR & 7);  // 1 & 63 gives 3.358 // shuttle 31 and 15: 3.250 //ece 31 and 1: 3.252
       /// These need to be masked   // shuttle : 1023 63 and 15: 3.254 // ece 1023, 31 and 1 : 3.254
       ////  63 and 7  and PC  -2 -4 -6 // 63 and 1 gives 3.256  // 63 and 7 s: // 31 and 7 : 3.243 best !
       for(int i = 0; i < NUMTAGTABLES; i++)
       {
            indexTagPred[i] = indexTagPred[i] & index_mask;
           
       }
       
        // get two predictions prime and alt (alternate)
       primePred = -1;
       altPred = -1;
       primeBank = NUMTAGTABLES;
       altBank = NUMTAGTABLES;
       
       /// See if any tag matches
       // T0 with longest history so if hit that awesome
        
       for(int iterator = 0; iterator < NUMTAGTABLES; iterator++)
       {
           
           
            if(tagPred[iterator][indexTagPred[iterator]].tag == tag[iterator])
            {
                primeBank = iterator;
                break;
            }  
       }      
            for(int iterator = primeBank + 1; iterator < NUMTAGTABLES; iterator++)
            {
                if(tagPred[iterator][indexTagPred[iterator]].tag == tag[iterator])
                {
                    altBank = iterator;
                    break;
                }  
            }    
            
       
       
    if(primeBank < NUMTAGTABLES)
    {        
       if(altBank == NUMTAGTABLES)
       {
           altPred = basePrediction;
       }
       else
       {
           if(tagPred[altBank][indexTagPred[altBank]].ctr >= TAGPRED_CTR_MAX/2)
                altPred = TAKEN;
            else 
                altPred = NOT_TAKEN;
       }
        
        if((tagPred[primeBank][indexTagPred[primeBank]].ctr  != 3) ||(tagPred[primeBank][indexTagPred[primeBank]].ctr != 4 ) || (tagPred[primeBank][indexTagPred[primeBank]].usefulBits != 0) || (altBetterCount < 8))
        {
            if(tagPred[primeBank][indexTagPred[primeBank]].ctr >= TAGPRED_CTR_MAX/2)
                primePred = TAKEN;
            else 
                primePred = NOT_TAKEN;
            return primePred;
        }
        else
        {
            return altPred;
        }
    }
    else
    {
        altPred = basePrediction;
        return altPred;
    }
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
 
void  PREDICTOR::UpdatePredictor(UINT32 PC, bool resolveDir, bool predDir, UINT32 branchTarget){
 
    bool strong_old_present = false;
    bool new_entry = 0;    
    if (primeBank < NUMTAGTABLES)
    {
        // As per update policy
        // 1st update the useful counter
        if ((predDir != altPred))
        {
	    
	    if (predDir == resolveDir)
	    {

		tagPred[primeBank][indexTagPred[primeBank]].usefulBits = SatIncrement(tagPred[primeBank][indexTagPred[primeBank]].usefulBits, BIMODAL_CTR_MAX);

	    }
	    else
	    {
		tagPred[primeBank][indexTagPred[primeBank]].usefulBits = SatDecrement(tagPred[primeBank][indexTagPred[primeBank]].usefulBits);
	    }

	}    
	 // 2nd update the counters which provided the prediction  
        if(resolveDir)
        {
            tagPred[primeBank][indexTagPred[primeBank]].ctr = SatIncrement(tagPred[primeBank][indexTagPred[primeBank]].ctr, TAGPRED_CTR_MAX);
        }
        else
        {
            tagPred[primeBank][indexTagPred[primeBank]].ctr = SatDecrement(tagPred[primeBank][indexTagPred[primeBank]].ctr);
        }
    }
    else
    {
        UINT32 bimodalIndex   = (PC) % (numBimodalEntries);
        if(resolveDir)
        {
            bimodal[bimodalIndex] = SatIncrement(bimodal[bimodalIndex], BIMODAL_CTR_MAX);
        }
        else
        {
            bimodal[bimodalIndex] = SatDecrement(bimodal[bimodalIndex]);
        }
    }
    // Check if the current Entry which gave the prediction is a newly allocated entry or not.
	if (primeBank < NUMTAGTABLES)
	{
	    
	    if((tagPred[primeBank][indexTagPred[primeBank]].usefulBits == 0) &&((tagPred[primeBank][indexTagPred[primeBank]].ctr  == 3) || (tagPred[primeBank][indexTagPred[primeBank]].ctr  == 4)))
            {
                new_entry = true;

		if (primePred != altPred)
		  {
		    if (altPred == resolveDir)
		      {
// Alternate prediction more useful is a counter to be of 4 bits
			if (altBetterCount < 15)
			{  
                            altBetterCount++;
                        }    
		      }

		    else if (altBetterCount > 0)
		    {
                        altBetterCount--;
                    }
                }
	    }
	}


// Proceeding to allocation of the entry.
    if((!new_entry) || (new_entry && (primePred != resolveDir)))
    {    
	if (((predDir != resolveDir) & (primeBank > 0)))
	  {
                    
	    for (int i = 0; i < primeBank; i++)
	      {
		if (tagPred[i][indexTagPred[i]].usefulBits == 0);
                strong_old_present = true;

	      }
// If no entry useful than decrease useful bits of all entries but do not allocate
	    if (strong_old_present == false)
	    {
		for (int i = primeBank - 1; i >= 0; i--)
		{
		    tagPred[i][indexTagPred[i]].usefulBits--;
                }
            }
	    else
	      {

                srand(time(NULL));
                int randNo = rand() % 100;
                int count = 0;
                int bank_store[NUMTAGTABLES - 1] = {-1, -1, -1};
                int matchBank = 0;
                for (int i = 0; i < primeBank; i++)
                {
                    if (tagPred[i][indexTagPred[i]].usefulBits == 0)
                    {
                        count++;
                        bank_store[i] = i;
                    }
                }  
                if(count == 1)
                {
                    matchBank = bank_store[0];
                }
                else if(count > 1)
                {
                     if(randNo > 33 && randNo <= 99)
                    {
                        matchBank = bank_store[(count-1)];
                    }
                    else
                    {
                        matchBank = bank_store[(count-2)];
                    }   
                }


		for (int i = matchBank; i > -1; i--)
		{
		    if ((tagPred[i][indexTagPred[i]].usefulBits == 0))
		      {
                        if(resolveDir)
                        {    
                            tagPred[i][indexTagPred[i]].ctr = 4;
                        }
                        else
                        {
                            tagPred[i][indexTagPred[i]].ctr = 3;
                        }    
                            tagPred[i][indexTagPred[i]].tag = tag[i];
                            tagPred[i][indexTagPred[i]].usefulBits = 0;
			break; // Only 1 entry allocated
		     }
                }

	    }

	}
    }    


// Periodic Useful bit Reset Logic ( Important so as to optimize compared to PPM paper)
	clock++;
        //for every 256 K instruction 1st MSB than LSB
	if(clock == (256*1024))
        {
            // reset clock
            clock = 0;
            if(clock_flip == 1)
            {
                // this is the 1st time
                clock_flip = 0;
            }
            else
            {
                clock_flip = 1;
            }
	    if(clock_flip == 1)
            {// MSB turn
                for (int jj = 0; jj < NUMTAGTABLES; jj++)
                {    
                    for (UINT32 ii = 0; ii < numTagPredEntries; ii++)
                    {
                        tagPred[jj][ii].usefulBits = tagPred[jj][ii].usefulBits & 1;
                    }
                }
            }    
            else
            {// LSB turn
                for (int jj = 0; jj < NUMTAGTABLES; jj++)
                   {    
                       for (UINT32 ii = 0; ii < numTagPredEntries; ii++)
                       {
                           tagPred[jj][ii].usefulBits = tagPred[jj][ii].usefulBits & 2;
                       }
                   }

	    }

	}

	
  // update the GHR
  GHR = (GHR << 1);

  if(resolveDir == TAKEN){
    GHR.set(0,1); 
  }

    for (int i = 0; i < NUMTAGTABLES; i++)
    {
                
            indexComp[i].updateCompHist(GHR);
            tagComp[0][i].updateCompHist(GHR);
            tagComp[1][i].updateCompHist(GHR);
            
    }
  // PHR update is simple, jus take the last bit ??
    // Always Limited to 16 bits as per paper.
    PHR = (PHR << 1); 
    if(PC & 1)
    {
        PHR = PHR + 1;
    }
    PHR = (PHR & ((1 << 16) - 1));
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void    PREDICTOR::TrackOtherInst(UINT32 PC, OpType opType, UINT32 branchTarget){

  // This function is called for instructions which are not
  // conditional branches, just in case someone decides to design
  // a predictor that uses information from such instructions.
  // We expect most contestants to leave this function untouched.

  return;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
