/* IRT scale score transformation

This file, which is part of Equating Recipes, is free softrware.
You can distribute it and/or modify it under the terms of the
GNU Lesser General Public License, version 3, as published by 
the Free Software Foundation.

This file is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License, version 3, for more details.

You should have received a copy of the GNU Lesser General Public 
License, version 3, in a ReadMe file distributed along with this
file.  If not, see <http://www.gnu.org/licenses/>   

Copyright 2009 
Center for Advanced Studies in Measurement and Assessment (CASMA)
University of Iowa

   Note: All pointers are 0-offset unless otherwise indicated.
  
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "IRTst.h" 


static struct IRTstControl *ContHandle;      /* global setting to control the method*/
static struct CommonItemSpec *ContComItem;
static enum symmetry ContSym;
static enum OnOff ContFuncStd; 

/***********************************************************************************/           

void ScaleTransform(const char *ItemOutF, const char *DistOutF, double slope,
		double intercept, struct ItemSpec *NewItem, struct IRTstControl *Handle)
/*--------------------------------------------------------------------------------------
  Functionality:
    Use the values of the slope and intercept of the new-to-old transformation
    to convert both the parameter estimates of the new form items and
    the ability points for the new group distribution.
    
    Input:
    
      ItemOutF  The name of a file in which the transformed output for items on
                the new form is saved; The format of output is, in essence, the
                same as that of input, so the output file can be read by the function
                ItemInfoRead without any syntax error.
      DistOutF  The name of a file in which the transformed ability points for the
                new group distribution are saved along with the original weights;
                As with ItemOutF, the output is saved using the format of input for
                the ability distribution, so the file DistOutF can be read by the
                function ThetaInfoRead without problems.
      slope     The value of A (slope) for a chosen linking method
      intercept The value of B (intercept) for a chosen linking method
      NewItem   A pointer to an array of the ItemSpec structure, which is for the new form
      Handle    A pointer to a variable of the IRTstControl structure
    
    Output:
      ItemOutF  The transformed output file for the new form
      DistOutF  The transformed output file for the new group's ability distribution

      * Note: The original input remains intact for both the items and distribution.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
--------------------------------------------------------------------------------------*/   
{
	int i, k;
	FILE *outf;

	/* open the output file for the transformed item parameters */
	outf=fopen(ItemOutF, "w");
	if(!outf) 
		runerror("can't open output file for the transformed item parameters\n");

	/* save the transformed information for the new form items */
	fprintf(outf, "%d\n", Handle->NewItemNum);

	for (i = 1; i <= Handle->NewItemNum; i++) {
		fprintf(outf, "%3d", NewItem[i].ItemID); /* write Item ID */	
		switch(NewItem[i].model)
		{
			case l3:
				fprintf(outf, "  L3 2 0 1"); /* write model, category number and score functions */
				fprintf(outf, " %g\n", NewItem[i].ScaleConst); /* write scaling constant */
				fprintf(outf, " %9.5f", NewItem[i].a[2]/slope);    /* write t-a (tranformed-a) */
				fprintf(outf, " %9.5f", NewItem[i].b[2]*slope + intercept); /* write t-b */
				fprintf(outf, " %9.5f", NewItem[i].c[2]);         /* write original c */        
				break;
			case gr:
				fprintf(outf, "  GR %d", NewItem[i].CatNum); /* write model and category number */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					fprintf(outf, " %g", NewItem[i].ScoreFunc[k]); /* write scoring function */
				}
				fprintf(outf, " %g\n", NewItem[i].ScaleConst); /* write scaling constant */
				fprintf(outf, " %9.5f", NewItem[i].a[2]/slope); /* write t-a */
				for(k=2; k <= NewItem[i].CatNum; k++) {
					if(k%10 == 1) fscanf(outf, "\n");
					fprintf(outf, " %9.5f", NewItem[i].b[k]*slope + intercept); /* write t-b */
				}
				break;
			case pc:
				fprintf(outf, "  PC %d", NewItem[i].CatNum); /* write model and category number */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					fprintf(outf, " %g", NewItem[i].ScoreFunc[k]); /* write scoring function */
				}
				fprintf(outf, " %g\n", NewItem[i].ScaleConst); /* write scaling constant */
				fprintf(outf, " %9.5f", NewItem[i].a[2]/slope); /* write t-a */
				fprintf(outf, " %9.5f", NewItem[i].b[0]*slope + intercept); /* write t-b */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					if ((k+2)%10 == 1) fprintf(outf, "\n");
					fprintf(outf, " %9.5f", NewItem[i].d[k]*slope); /* write t-d */
				}
				break;
			case nr:
				fprintf(outf, "  NR %d", NewItem[i].CatNum); /* write model and category number */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					fprintf(outf, " %g", NewItem[i].ScoreFunc[k]); /* write scoring function */
				}
				fprintf(outf, " 1.0\n"); /* write scaling constant */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					if (k > 1 && k%10 == 1) fprintf(outf, "\n");
					fprintf(outf, " %9.5f", NewItem[i].a[k]/slope); /* write t-a */
				}
				for(k=1; k <= NewItem[i].CatNum; k++) {
					if (k%10 == 1) fprintf(outf, "\n");
					/* write t-c */
					fprintf(outf, " %9.5f", NewItem[i].c[k] - (intercept/slope)*NewItem[i].a[k]);
				}
				break;
			default:
				break;
		} /* end of switch */
		fprintf(outf, "\n");
	} /* end of for loop */
	fclose(outf);


	/* open the output file for the transformed distribution for the new group */
	outf=fopen(DistOutF, "w");
	if(!outf) 
		runerror("can't open output file for the transformed"
                 " distribution for the new group\n");

	/* save the transformed information for the new group */
	fprintf(outf, "%d\n", Handle->NewThetaNum);
	
	for (i = 1; i <= Handle->NewThetaNum; i++) {
		fprintf(outf, "%E", Handle->NewThetaValues[i]*slope + intercept);
		fprintf(outf, "   %E\n", Handle->NewThetaWeights[i]);
	}
	fclose(outf);

	return;
}

/***********************************************************************************/
 
struct ItemSpec *ItemInfoRead(FILE *inf, const char *oldOrnew, struct IRTstControl *Handle)
/*-----------------------------------------------------------------------------------------
  Functionality:
    Read pieces of information for items on the new or old form from an input
    file, depending on the oldOrnew string, which must be "new" or "old."
    Assume that the input file is opened but is not read in yet at all.
    To read the information, a pointer to the ItemSpec structure is first created
    using the function malloc, and then some of the structure members,
    if they are pointers, are initialized with designated memory addresses using
    the function malloc. In addition, some members of an IRTstControl structure
    referenced by the pointer Handle are initialized to their respective values
    that have been read.
    
    Input:
      inf       A file pointer; the file must be opened for input.
      oldOrnew  A string indicating that the items to be read are on the old or new form.
                The string must be either "new" or "old"
      Handle    A pointer to a variable of the IRTstControl structure, which is defined
                in the header file IRTst.h

      The input file must be prepared according to the following format:
    
      Line 1    Number of items on the particular form (new or old)
      Line 2+   A record for each item is typed. The record can range over multiple
                lines. So, fscanf() is used, instead of using fgets().
                The general format of the record is
    	      
    	        <ID>  <Model> <CatNum> <ScoreFunc> <ScaleConst> <ItemParsList>
    	        --------------------------------------------------------------
    	        ID            Item ID
    	        Model         L3 (three-parameter logistic model)
    	        (Uppercase)   GR (logistic graded response model)
    	                      PC (generalized partial credit model)
    	                      NR (nominal response model)
    	        CatNum        Number of categories for the item
    	        ScoreFunc     A list of scores given to CatNum categories
    	        ScaleConst    Scaling Constant
    	                      In the case of the NR model, it should equal 1;
    	                      if other values are typed for the NR model,
    	                      they are ignored.
    	                    
    	        ItemParsList  A list of item parameters for the item
    	        --------------------------------------------------------------
    	      
    	        Specifically, for example, the item record is typed differently
    	        according to the model used as follows (assume that item j has
    	        K categories). Note that curly braces are used to group score
    	        functions and item parameters, but they must not be typed when
    	        preparing the input file: 
    
                1. The three-parameter logistic model (L3):
              
     	           ID L3 2 {0 1} 1.7 {a b c}
     	      
     	        2. The graded response model (GR):
     	      
                   ID GR 4 {0 1 2 3} 1.7 {a b2 b3 b4}
                 
                   * Note: b1 parameter does not exist under the GR model
              
                3. The generalized partial credit model (PC):
              
                   ID PC 4 {1 2 3 4} 1.7 {a b d1 d2 d3 d4}
                 
                   * Note: Location paramter b and category parameters d1 - d4 are
                           used, instead of using the item-step parameter bd,
                           where bd = b - d.
                4. The nominal response model (NR):
              
                   ID NR 3 {0 1 2} 1.0 {a1 a2 a3 b1 b2 b3}
                 
                   * Note: Assume that category responses are ordered.
    
    Output:
      Return a pointer to the ItemSpec structure.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
-----------------------------------------------------------------------------------------*/
    
{
	char buff[3];
	int i, k, ItemNum;
	struct ItemSpec *Item;

	fscanf(inf, "%d",   &ItemNum); /* reading the number of items */
	if (strcmp(oldOrnew, "old") == 0) Handle->OldItemNum = ItemNum;
	else                              Handle->NewItemNum = ItemNum;

	Item = (struct ItemSpec *) malloc((ItemNum+1)*sizeof(struct ItemSpec));
	if (!Item) 
		runerror("memory allocation failure 1 in ItemInfoRead()\n");

	for (i = 1; i <= ItemNum; i++) {
		fscanf(inf, "%d", &(Item[i].ItemID)); /* read ID */
		
		fscanf(inf, "%s", buff);               /* read model */
		if (strcmp(buff, "L3")==0 || strcmp(buff, "l3")==0)      Item[i].model = l3;
		else if (strcmp(buff, "GR")==0 || strcmp(buff, "gr")==0) Item[i].model = gr;
		else if (strcmp(buff, "PC")==0 || strcmp(buff, "pc")==0) Item[i].model = pc;
		else if (strcmp(buff, "NR")==0 || strcmp(buff, "nr")==0) Item[i].model = nr;
        else runerror("Invalid value for type of IRT model");

	
		fscanf(inf, "%d", &(Item[i].CatNum)); /* read CatNum */
                
		Item[i].ScoreFunc = (double *) malloc((Item[i].CatNum+1) * sizeof(double));
		if (!Item[i].ScoreFunc) 
			runerror("memory allocation failure 2 in ItemInfoRead()\n");
                
		for (k = 1; k <= Item[i].CatNum; k++) { /* read scoring function */
			fscanf(inf, "%lf", &(Item[i].ScoreFunc[k]));
		}
		
                
		Item[i].a = (double *) malloc((Item[i].CatNum+1)*sizeof(double));
		if (!Item[i].a) 
			runerror("memory allocation failure 3 in ItemInfoRead()\n");

		Item[i].b = (double *) malloc((Item[i].CatNum+1)*sizeof(double));
		if (!Item[i].b) 
			runerror("memory allocation failure 4 in ItemInfoRead()\n");

		Item[i].c = (double *) malloc((Item[i].CatNum+1)*sizeof(double));
		if (!Item[i].c) 
			runerror("memory allocation failure 5 in ItemInfoRead()\n");

		Item[i].d = (double *) malloc((Item[i].CatNum+1)*sizeof(double));
		if (!Item[i].d) 
			runerror("memory allocation failure 6 in ItemInfoRead()\n");

                		
		switch(Item[i].model)
		{
			case l3:
				fscanf(inf, "%lf", &(Item[i].ScaleConst));   /* read D */
				fscanf(inf, "%lf", &(Item[i].a[2]));         /* read a parameter */
				fscanf(inf, "%lf", &(Item[i].b[2]));         /* read b parameter */
				fscanf(inf, "%lf", &(Item[i].c[2]));         /* read c parameter */        
				break;
			case gr:
				fscanf(inf, "%lf", &(Item[i].ScaleConst));   /* read D */
				fscanf(inf, "%lf", &(Item[i].a[2]));         /* read a parameter */
				for(k=2; k <= Item[i].CatNum; k++) {
					fscanf(inf, "%lf", &(Item[i].b[k])); /* read b parameters */
				}
				break;
			case pc:
				fscanf(inf, "%lf", &(Item[i].ScaleConst));   /* read D */
				fscanf(inf, "%lf", &(Item[i].a[2]));         /* read a parameter */
				fscanf(inf, "%lf", &(Item[i].b[0]));         /* read b parameter */
				for(k=1; k <= Item[i].CatNum; k++) {
					fscanf(inf, "%lf", &(Item[i].d[k])); /* read d parameters */
					Item[i].b[k] = Item[i].b[0] - Item[i].d[k];
				}
				break;
			case nr:
				fscanf(inf, "%lf", &(Item[i].ScaleConst));   /* read D = 1.0 */
				if (Item[i].ScaleConst != 1.0) Item[i].ScaleConst = 1.0;
				for(k=1; k <= Item[i].CatNum; k++) {
					fscanf(inf, "%lf", &(Item[i].a[k])); /* read a parameters */
				}
				for(k=1; k <= Item[i].CatNum; k++) {
					fscanf(inf, "%lf", &(Item[i].c[k])); /* read c parameters */
				}
				break;
			default:
				break;
		} /* end of switch */
	} /* end of for loop */      
	return Item;
}

/***********************************************************************************/
 
struct CommonItemSpec *ItemPairsRead(FILE *inf, struct ItemSpec *NewItem,
		struct ItemSpec *OldItem, struct IRTstControl *Handle)
/*------------------------------------------------------------------------------
  Functionality:
    Read pairs of new and old items, which are common items, from an input file.
    Assume that the input file is opened but is not read in yet at all.
    To read the information, a pointer to the ItemSpec structure is first created
    using the function malloc. While reading the pairs information, an array that
    is referenced by the pointer is set up with the two ItemSpec "arrays," NewItem
    and OldItem. The common information between the new and old items, such as
    ScaleConst, CatNum, and model, is set with the new item information.
    In addition, the ComItemNum member of an IRTstControl structure referenced
    by the pointer Handle is initialized to the number of common items.

   Input:
      inf       A file pointer; the file must be opened for input.
      NewItem   A pointer to the ItemSpec structure; the pointer points to the
                memory address of an array of the ItemSpec structure for items on
                the new form. The ItemSpec structure is defined in the header file
                IRTst.h
      OldItem   A pointer to the ItemSpec structure, which is used for items on
                the old form
      Handle    A pointer to a variable of the IRTstControl structure, which is defined
                in the header file IRTst.h
    
      The input file must be prepared according to the following format:
    
      Line 1    Number of the common items, which are used to estimate scaling
                constants, A and B.
      Line 2+   Pairs of the new and old item ID's.
                Each pair consists of two integers, a new item ID and
                a matching old item ID.

                Suppose, for example, that there are two common items and the
                first and second items on the new form match with the third and
                fourth items on the old form. Pairs of items, then, must be
                designated as follows, with the first line being for the number of
                common items:
         
                2
                1 3
                2 4

    Output:
      Return a pointer to the CommonItemSpec structure.              
      While reading the pairs information, the ComItem array is set with the two 
      ItemSpec arrays, NewItem and OldItem.
      The common information between the new and old items, such as ScaleConst,
      CatNum, and model, is set with the new item information.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int i, k, l, it1, it2, ItemNum, id_finding_error, Nid, Oid;
	struct CommonItemSpec *ComItem;

	fscanf(inf, "%d",   &ItemNum); /* reading the number of common items */
	if (ItemNum <= Handle->NewItemNum && ItemNum <= Handle->OldItemNum)
		Handle->ComItemNum = ItemNum;
	else 
		runerror("number of common items must not exceed that of new or old items");


	ComItem = (struct CommonItemSpec *) malloc((ItemNum+1)*sizeof(struct CommonItemSpec));
	if (!ComItem) 
		runerror("memory allocation failure 1 in ItemPairsRead()");

        
	for (i = 1; i <= ItemNum; i++) {
		fscanf(inf, "%d", &it1); /* read new item ID */
		fscanf(inf, "%d", &it2); /* read old item ID */
		
		ComItem[i].NewID = it1;
		ComItem[i].OldID = it2;

       /* finding memory indices for it1 and it2 */
		id_finding_error = 0;
		for(l = 1; l <= Handle->NewItemNum; l++) {
			if (NewItem[l].ItemID == it1) {
				Nid = l;
				id_finding_error ++;
			}
		}
		if (id_finding_error != 1) {
			fprintf(stderr, "matching error in common item %d: new item's ID cannot be found or is used more than once\n", i);
			runerror("Check the input data for common items\n");
		}

		id_finding_error = 0;
		for(l = 1; l <= Handle->OldItemNum; l++) {
			if (OldItem[l].ItemID == it2) {
				Oid = l;
				id_finding_error ++;
			}
		}
		if (id_finding_error != 1) {
			fprintf(stderr, "matching error in common item %d: old item's ID cannot be found or is used more than once\n", i);
			runerror("Check the input data for common items\n");
		}

		/* minimum checking for the item matching */
		if (NewItem[Nid].model == OldItem[Oid].model)
			ComItem[i].model = NewItem[Nid].model;
		else {
			fprintf(stderr, "matching error in IRT model between new item %d and old item %d\n", it1, it2);
			runerror("Check the input data for new and old items\n");
		}

		if (NewItem[Nid].CatNum == OldItem[Oid].CatNum)
			ComItem[i].CatNum = NewItem[Nid].CatNum;
		else {
			fprintf(stderr, "matching error in category number between new item %d and old item %d\n", it1, it2);
			runerror("Check the input data for new and old items\n");
		}
		
		if (NewItem[Nid].ScaleConst == OldItem[Oid].ScaleConst)
			ComItem[i].ScaleConst = NewItem[Nid].ScaleConst;
		else {
			fprintf(stderr, "matching error in scoring function between new item %d and old item %d\n", it1, it2);
			runerror("Check the input data for new and old items\n");
		}

		/* memory allocation for scoring function */                
		ComItem[i].ScoreFunc = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].ScoreFunc) 
			runerror("memory allocation failure 2 in ItemPairsRead()");

                
		for (k = 1; k <= ComItem[i].CatNum; k++) { /* assigning scoring function */
			ComItem[i].ScoreFunc[k] = NewItem[Nid].ScoreFunc[k];
		}
		
		/* memory allocation for item parameters from new and old forms */
		ComItem[i].Na = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].Na) 
			runerror("memory allocation failure 3 in ItemPairsRead()");

		ComItem[i].Nb = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].Nb) 
			runerror("memory allocation failure 4 in ItemPairsRead()");

		ComItem[i].Nc = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].Nc) 
			runerror("memory allocation failure 5 in ItemPairsRead()");

		ComItem[i].Nd = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].Nd) 
			runerror("memory allocation failure 6 in ItemPairsRead()");

                
		ComItem[i].Oa = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].Oa) 
			runerror("memory allocation failure 7 in ItemPairsRead()");

		ComItem[i].Ob = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].Ob) 
			runerror("memory allocation failure 8 in ItemPairsRead()");

		ComItem[i].Oc = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].Oc) 
			runerror("memory allocation failure 9 in ItemPairsRead()");

		ComItem[i].Od = (double *) malloc((ComItem[i].CatNum+1)*sizeof(double));
		if (!ComItem[i].Od) 
			runerror("memory allocation failure 10 in ItemPairsRead()");

                		
		switch(ComItem[i].model)
		{
			case l3:
				ComItem[i].Na[2] = NewItem[Nid].a[2];
				ComItem[i].Nb[2] = NewItem[Nid].b[2];
				ComItem[i].Nc[2] = NewItem[Nid].c[2];
				ComItem[i].Oa[2] = OldItem[Oid].a[2];
				ComItem[i].Ob[2] = OldItem[Oid].b[2];
				ComItem[i].Oc[2] = OldItem[Oid].c[2];
				break;
			case gr:
				ComItem[i].Na[2] = NewItem[Nid].a[2];
				ComItem[i].Oa[2] = OldItem[Oid].a[2];
				for (k=2; k <= ComItem[i].CatNum; k++) {
					ComItem[i].Nb[k] = NewItem[Nid].b[k];
					ComItem[i].Ob[k] = OldItem[Oid].b[k];
				}
				break;
			case pc:
				ComItem[i].Na[2] = NewItem[Nid].a[2];
				ComItem[i].Nb[0] = NewItem[Nid].b[0];
				ComItem[i].Oa[2] = OldItem[Oid].a[2];
				ComItem[i].Ob[0] = OldItem[Oid].b[0];
				for (k=1; k <= ComItem[i].CatNum; k++) {
					ComItem[i].Nb[k] = NewItem[Nid].b[k];
					ComItem[i].Nd[k] = NewItem[Nid].d[k];
					ComItem[i].Ob[k] = OldItem[Oid].b[k];
					ComItem[i].Od[k] = OldItem[Oid].d[k];
				}
				break;
			case nr:
				for (k=1; k <= ComItem[i].CatNum; k++) {
					ComItem[i].Na[k] = NewItem[Nid].a[k];
					ComItem[i].Nc[k] = NewItem[Nid].c[k];
					ComItem[i].Oa[k] = OldItem[Oid].a[k];
					ComItem[i].Oc[k] = OldItem[Oid].c[k];
				}
				break;
			default:
				break;
		} /* end of switch */
	} /* end of for loop */      
	return ComItem;
}

/***********************************************************************************/

double ProbOld(struct CommonItemSpec *Item, int CatID, double theta, enum OnOff original,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Calculate the value of the characteristic curve on the old scale by model.

 For details about partial derivatives, refer to the following report:

   Kim, S., & Kolen, M.J. (2005). Methods for obtaining a common scale under
      unidimensional IRT models: A technical review and further extensions
      (Iowa Testing Programs Occasional Paper, No. 52). The University of
      Iowa.

  Input:
	Item: One item of struct CommonItemSpec type 
	CatID: Response category ID in question (1 through CatNum)
	theta: ability value
	original: on or off
		if on, then old scale's item parameters are used with
		S = 1.0 and I = 0.0. In this case, S and I are over-argumented.
		Otherwise, transformed item parameters (from new scale) are
		used through S and I.
	S: slope of the linear tranformation
	I: intercept of the linear transformation

  Output:
	Return either original (with old parameters and S = 1.0 and I = 0.0)
	or transformed (with new parameters and S and I) probability
	on the old scale

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double prob;

	switch(Item->model)
	{
		case l3:
			if (original == on)
				prob = Prob3PL(CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob[2], Item->Oc[2], "old", 1, 0);
			else
				prob = Prob3PL(CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb[2], Item->Nc[2], "old", S, I);
			break;
		case gr:
			if (original == on)
				prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob, "old", 1, 0);
			else
				prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb, "old", S, I);
			break;
		case pc:
			if (original == on)
				prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob, "old", 1, 0);
			else
				prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb, "old", S, I);
			break;
		case nr:
			if (original == on)
				prob = ProbNRM(Item->CatNum, CatID, theta,
					Item->Oa, Item->Oc, "old", 1, 0);
			else
				prob = ProbNRM(Item->CatNum, CatID, theta,
					Item->Na, Item->Nc, "old", S, I);
			break;
		default:
			break;
	}
	return prob;
}

/***********************************************************************************/

double ProbNew(struct CommonItemSpec *Item, int CatID, double theta, enum OnOff original,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Calculate the value of the characteristic curve on the new scale by model.

    Input:
	Item: One item of struct CommonItemSpec type 
	CatID: Response category ID in question (1 through CatNum)
	theta: ability value
	original: on or off
		if on, then new scale's item parameters are used with
		S = 1.0 and I = 0.0. In this case, S and I are over-argumented.
		Otherwise, transformed item parameters (from old scale) are
		used through S and I.
	S: slope of the linear tranformation
	I: intercept of the linear transformation

   Output:
	Return either original (with new parameters and S = 1.0 and I = 0.0)
	or transformed (with old parameters and S and I) probability
	on the new scale

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double prob;

	switch(Item->model)
	{
		case l3:
			if (original == on)
				prob = Prob3PL(CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb[2], Item->Nc[2], "new", 1, 0);
			else
				prob = Prob3PL(CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob[2], Item->Oc[2], "new", S, I);
			break;
		case gr:
			if (original == on)
				prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb, "new", 1, 0);
			else
				prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob, "new", S, I);
			break;
		case pc:
			if (original == on)
				prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb, "new", 1, 0);
			else
				prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
                        		Item->Oa[2], Item->Ob, "new", S, I);
			break;
		case nr:
			if (original == on)
				prob = ProbNRM(Item->CatNum, CatID, theta,
					Item->Na, Item->Nc, "new", 1, 0);
			else
				prob = ProbNRM(Item->CatNum, CatID, theta,
					Item->Oa, Item->Oc, "new", S, I);
			break;
		default:
			break;
	}
	return prob;
}

/***********************************************************************************/

double PdOldOverS(struct CommonItemSpec *Item, int CatID, double theta, double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    By model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to S.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double pd;

	switch(Item->model)
	{
		case l3:
			pd = Pd3PLOldOverS(CatID, theta, Item->ScaleConst,
				Item->Na[2], Item->Nb[2], Item->Nc[2], S, I);
			break;
		case gr:
			pd = PdLGROldOverS(Item->CatNum, CatID, theta, Item->ScaleConst,
                                Item->Na[2], Item->Nb, S, I);
			break;
		case pc:
			pd = PdGPCOldOverS(Item->CatNum, CatID, theta, Item->ScaleConst,
                                Item->Na[2], Item->Nb, S, I);
			break;
		case nr:
			pd = PdNRMOldOverS(Item->CatNum, CatID, theta,
				Item->Na, Item->Nc, S, I);
			break;
		default:
			break;
	}
	return pd;
}

/***********************************************************************************/

double PdOldOverI(struct CommonItemSpec *Item, int CatID, double theta, double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    By model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to I.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double pd;

	switch(Item->model)
	{
		case l3:
			pd = Pd3PLOldOverI(CatID, theta, Item->ScaleConst,
				Item->Na[2], Item->Nb[2], Item->Nc[2], S, I);
			break;
		case gr:
			pd = PdLGROldOverI(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Na[2], Item->Nb, S, I);
			break;
		case pc:
			pd = PdGPCOldOverI(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Na[2], Item->Nb, S, I);
			break;
		case nr:
			pd = PdNRMOldOverI(Item->CatNum, CatID, theta,
				Item->Na, Item->Nc, S, I);
			break;
		default:
			break;
	}
	return pd;
}


/***********************************************************************************/

double PdNewOverS(struct CommonItemSpec *Item, int CatID, double theta, double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    By model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to S.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double pd;
   
	switch(Item->model)
	{
		case l3:
			pd = Pd3PLNewOverS(CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob[2], Item->Oc[2], S, I);
			break;
		case gr:
			pd = PdLGRNewOverS(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob, S, I);
			break;
		case pc:
			pd = PdGPCNewOverS(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob, S, I);
			break;
		case nr:
			pd = PdNRMNewOverS(Item->CatNum, CatID, theta,
				Item->Oa, Item->Oc, S, I);
			break;
		default:
			break;
	}
	return pd;
}

/***********************************************************************************/

double PdNewOverI(struct CommonItemSpec *Item, int CatID, double theta, double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    By model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to I.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double pd;

	switch(Item->model)
	{
		case l3:
			pd = Pd3PLNewOverI(CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob[2], Item->Oc[2], S, I);
			break;
		case gr:
			pd = PdLGRNewOverI(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob, S, I);
			break;
		case pc:
			pd = PdGPCNewOverI(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob, S, I);
			break;
		case nr:
			pd = PdNRMNewOverI(Item->CatNum, CatID, theta,
				Item->Oa, Item->Oc, S, I);
			break;
		default:
			break;
	}
	return pd;
}

/***********************************************************************************/

double Prob3PL(int CatID, double theta, double D, double a, double b, double c,
		char scale[], double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Calculate the value of the characteristic curve for the 3PL model.
    Input: 
	CatID: response category, 1 or 2 (1 for incorrect; 2 for correct)
	theta: ability value
	D: scaling constant (typically 1.7)
	a, b, c: item parameters
	scale: "old" or "new" ability scale, on which item and ability parameter
		estimates are placed on.
	S: slope of the linear tranformation
	I: intercept of the linear transformation

   Output:
	Return either original or transformed probability
	For the original probability, S = 1 and I = 0 with parameters on the
	reference scale.
	For the transformed probability, S and I with parameters on the
	transformed scale.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as, bs, cs; /* parameters from new-to-old transformation */
	double ar, br, cr; /* parameters from old-to-new transformation */
	double devs, devr;

	if (strcmp(scale, "old") == 0) {
		/* new-to-old scale transformation */
		as = a/S;
		bs = S*b + I;
		cs = c;
		devs = D*as*(theta-bs);
   		return ( cs*(CatID-1) + (1.0-cs)*exp(-devs*(2-CatID))/(1.0+exp(-devs)) );
	}
	else {
		/* old-to-new scale transformation */
		ar = S*a;
		br = (b-I)/S;
		cr = c;
		devr = D*ar*(theta - br);
		return ( cr*(CatID-1) + (1.0-cr)*exp(-devr*(2-CatID))/(1.0+exp(-devr)) );
	}
}

/***********************************************************************************/

double Pd3PLOldOverS(int CatID, double theta, double D, double na, double nb, double nc,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to S.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as, cs;
	double ps, qs; /* Probability with transformed item parameters */
	double ps_over_S;

	/* Scale Transformation */
	as = na / S;
	cs = nc;

	ps = Prob3PL(2, theta, D, na, nb, nc, "old", S, I);
	qs = 1.0 - ps;
	ps_over_S = -D*as*( (theta - I)/S )*(ps - cs)*qs / (1.0 - cs);
   
	if (CatID == 1) return -ps_over_S;
	else if (CatID == 2) return ps_over_S;
}

/***********************************************************************************/

double Pd3PLOldOverI(int CatID, double theta, double D, double na, double nb, double nc,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to I.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as, cs;
	double ps, qs; /* Probability with transformed item parameters */
	double ps_over_I;

	/* Scale Transformation */
	as = na / S;
	cs = nc;

	ps = Prob3PL(2, theta, D, na, nb, nc, "old", S, I);
	qs = 1.0 - ps;
	ps_over_I = -D*as*(ps - cs)*qs / (1.0 - cs);
   
	if (CatID == 1) return -ps_over_I;
	else if (CatID == 2) return ps_over_I;
}

/***********************************************************************************/

double Pd3PLNewOverS(int CatID, double theta, double D, double oa, double ob, double oc,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to S.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double cr;
	double pr, qr; /* Probability with transformed item parameter estimates */
	double pr_over_S;

	/* Scale Transformation */
	cr = oc;

	pr = Prob3PL(2, theta, D, oa, ob, oc, "new", S, I);
	qr = 1.0 - pr;

	pr_over_S = D*oa*theta*(pr - cr)*qr / (1.0 - cr);
   
	if (CatID == 1) return -pr_over_S;
	else if (CatID == 2) return pr_over_S;
}

/***********************************************************************************/

double Pd3PLNewOverI(int CatID, double theta, double D, double oa, double ob, double oc,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to I.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
   double cr;
   double pr, qr; /* Probability with transformed item parameter estimates */
   double pr_over_I;

	/* Scale Transformation */
	cr = oc;

	pr = Prob3PL(2, theta, D, oa, ob, oc, "new", S, I);
	qr = 1.0 - pr;

	pr_over_I = D*oa*(pr - cr)*qr / (1.0 - cr);

	if (CatID == 1) return -pr_over_I;
	else if (CatID == 2) return pr_over_I;
}

/***********************************************************************************/

double CumProbLGR(int CatID, double theta, double D, double a, double b,
		char scale[], double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the logistic graded respnse model, calculate the cumulative category
    characteristic curve. Use the notation used in Kim and Kolen (2005).

    Input
	scale: "old" or "new"
    Output
	Return either original or transformed probability.
	For the original probability, S = 1 and I = 0 with parameters on the
	reference scale.
	For the transfomred probability, S and I with parameters on the
	transformed scale.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as, bs, ar, br;
	double cprob;

	if (CatID == 1) cprob = 1.0;
	else {
		if (strcmp(scale, "old") == 0) {
			/* new-to-old transformation */
			as = a/S;
			bs = S*b + I;
			cprob = 1.0/(1.0 + exp(-D*as*(theta - bs)));
		}
		else {
			/* old-to-new transformation */
			ar = S*a;
			br = (b - I)/S;
			cprob = 1.0/(1.0 + exp(-D*ar*(theta - br)));
		}
	}
	return cprob;
}

/***********************************************************************************/

double ProbLGR(int CatNum, int CatID, double theta, double D, double a, double b[],
               char scale[], double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the logistic graded response model, calculate the category
    characteristic curve. Use the notation used in Kim and Kolen (2005).

    Input
	CatNum: number of categories
	CatID: category response ID
	theta: ability value
	D: scaling constant
	a: discrimination parameter
	b[]: difficulty parameter array
	     b[2] for the first difficulty parameter
	     b[CatNum] for the last category
	scale: "old" or "new"
	S: slope of linear transformation
	I: intercept of linear transformation
    Output
	Return either original or transformed probability.
	For the original probability, S = 1 and I = 0 with parameters on the
	reference scale.
	For the transfomred probability, S and I with parameters on the
	transformed scale.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double pre_cp, pos_cp;
  
	if (strcmp(scale, "old") == 0) {
		pre_cp = CumProbLGR(CatID, theta, D, a, b[CatID], "old", S, I);
		if (CatID < CatNum)
			pos_cp = CumProbLGR(CatID+1, theta, D, a, b[CatID+1], "old", S, I);
		else
			pos_cp = 0.0;
	}
	else {
		pre_cp = CumProbLGR(CatID, theta, D, a, b[CatID], "new", S, I);
		if (CatID < CatNum) 
			pos_cp = CumProbLGR(CatID+1, theta, D, a, b[CatID+1], "new", S, I);
		else 
			pos_cp = 0.0;
	}
	return (pre_cp - pos_cp);
}

/***********************************************************************************/

double PdLGROldOverS(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the GR model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to S.
    
    Note:
	nb[2] for the first item-step parameter
	nb[CatNum] for the last item-step parameter

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as;
	double pre_cps, pos_cps;
	double pre_cps_over_S, pos_cps_over_S;
	double ps_over_S;

	as = na/S;

	if (CatID == 1) {
		pos_cps = CumProbLGR(CatID+1, theta, D, na, nb[CatID+1], "old", S, I);
		pos_cps_over_S = -D*as*((theta - I)/S)*pos_cps*(1.0 - pos_cps);
		ps_over_S = 0.0 - pos_cps_over_S;
	}
	else {
		if (CatID < CatNum) {
			pre_cps = CumProbLGR(CatID, theta, D, na, nb[CatID], "old", S, I);
			pos_cps = CumProbLGR(CatID+1, theta, D, na, nb[CatID+1], "old", S, I);
			pre_cps_over_S = -D*as*((theta - I)/S)*pre_cps*(1.0 - pre_cps);
			pos_cps_over_S = -D*as*((theta - I)/S)*pos_cps*(1.0 - pos_cps);
			ps_over_S = pre_cps_over_S - pos_cps_over_S;
		}
		else { /* CatId == CatNum */
			pre_cps = CumProbLGR(CatID, theta, D, na, nb[CatID], "old", S, I);
			pre_cps_over_S = -D*as*((theta - I)/S)*pre_cps*(1.0 - pre_cps);
			ps_over_S = pre_cps_over_S - 0.0;
		}
	}
	return ps_over_S;
}

/***********************************************************************************/

double PdLGROldOverI(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the GR model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to I.
    
    Note:
	nb[2] for the first item-step parameter
	nb[CatNum] for the last item-step parameter

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as;
	double pre_cps, pos_cps;
	double pre_cps_over_I, pos_cps_over_I;
	double ps_over_I;

	as = na / S;

	if (CatID == 1) {
		pos_cps = CumProbLGR(CatID+1, theta, D, na, nb[CatID+1], "old", S, I);
		pos_cps_over_I = -D*as*pos_cps*(1.0 - pos_cps);
		ps_over_I = 0.0 - pos_cps_over_I;
	}
	else {
		if (CatID < CatNum) {
			pre_cps = CumProbLGR(CatID, theta, D, na, nb[CatID], "old", S, I);
			pos_cps = CumProbLGR(CatID+1, theta, D, na, nb[CatID+1], "old", S, I);
			pre_cps_over_I = -D*as*pre_cps*(1.0 - pre_cps);
			pos_cps_over_I = -D*as*pos_cps*(1.0 - pos_cps);
			ps_over_I = pre_cps_over_I - pos_cps_over_I;
		}
		else { /* CatID == CatNum */
			pre_cps = CumProbLGR(CatID, theta, D, na, nb[CatID], "old", S, I);
			pre_cps_over_I = -D*as*pre_cps*(1.0 - pre_cps);
			ps_over_I = pre_cps_over_I - 0.0;
		}
	}
	return ps_over_I;
}

/***********************************************************************************/

double PdLGRNewOverS(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the GR model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to S.
    
    Note:
	ob[2] for the first item-step parameter
	ob[CatNum] for the last item-step parameter

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double pre_cpr, pos_cpr;
	double pre_cpr_over_S, pos_cpr_over_S;
	double pr_over_S;
  
	if (CatID == 1) {
		pos_cpr = CumProbLGR(CatID+1, theta, D, oa, ob[CatID+1], "new", S, I);
		pos_cpr_over_S = D*oa*theta*pos_cpr*(1.0 - pos_cpr);
		pr_over_S = 0.0 - pos_cpr_over_S;
	}
	else {
		if (CatID < CatNum) {
			pre_cpr = CumProbLGR(CatID, theta, D, oa, ob[CatID], "new", S, I);
			pos_cpr = CumProbLGR(CatID+1, theta, D, oa, ob[CatID+1], "new", S, I);
			pre_cpr_over_S = D*oa*theta*pre_cpr*(1.0 - pre_cpr);
			pos_cpr_over_S = D*oa*theta*pos_cpr*(1.0 - pos_cpr);
			pr_over_S = pre_cpr_over_S - pos_cpr_over_S;
		}
		else { /* CatID == CatNum */
			pre_cpr = CumProbLGR(CatID, theta, D, oa, ob[CatID], "new", S, I);
			pre_cpr_over_S = D*oa*theta*pre_cpr*(1.0 - pre_cpr);
			pr_over_S = pre_cpr_over_S - 0.0;
		}
	}
	return pr_over_S;
}

/***********************************************************************************/

double PdLGRNewOverI(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the GR model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to I.
    
    Note:
	ob[2] for the first item-step parameter
	ob[CatNum] for the last item-step parameter

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double pre_cpr, pos_cpr;
	double pre_cpr_over_I, pos_cpr_over_I;
	double pr_over_I;
  
	if (CatID == 1) {
		pos_cpr = CumProbLGR(CatID+1, theta, D, oa, ob[CatID+1], "new", S, I);
		pos_cpr_over_I = D*oa*pos_cpr*(1.0 - pos_cpr);
		pr_over_I = 0.0 - pos_cpr_over_I;
	}
	else {
		if (CatID < CatNum) {
			pre_cpr = CumProbLGR(CatID, theta, D, oa, ob[CatID], "new", S, I);
			pos_cpr = CumProbLGR(CatID+1, theta, D, oa, ob[CatID+1], "new", S, I);
			pre_cpr_over_I = D*oa*pre_cpr*(1.0 - pre_cpr);
			pos_cpr_over_I = D*oa*pos_cpr*(1.0 - pos_cpr);
			pr_over_I = pre_cpr_over_I - pos_cpr_over_I;
		}
		else { /* CatID == CatNum */
			pre_cpr = CumProbLGR(CatID, theta, D, oa, ob[CatID], "new", S, I);
			pre_cpr_over_I = D*oa*pre_cpr*(1.0 - pre_cpr);
			pr_over_I = pre_cpr_over_I - 0.0;
		}
	}
	return pr_over_I;
}

/***********************************************************************************/

double ProbGPC(int CatNum, int CatID, double theta, double D, double a, double b[],
		char scale[], double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the generalized partial credit model, calculate the category
    characteristic curve. Use the notation used in Kim and Kolen (2005).
    
    Note:
	is dealt with as a special case of the nominal response model.

    Input
	CatNum: number of categories
	CatID: category response ID
	theta: ability value
	D: scaling constant
	a: discrimination parameter
	b[]: difficulty parameter array
	     It is assumed that b[1] = 0, but is not used.
	     b[2] for the actual first difficulty (item-step) parameter
	     b[CatNum] for the last category
	scale: "old" or "new"
	S: slope of linear transformation
	I: intercept of linear transformation
    Output
	Return either original or transformed probability.
	For the original probability, S = 1 and I = 0 with parameters on the
	reference scale.
	For the transfomred probability, S and I with parameters on the
	transformed scale.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k, l;
	double vjs = 0.0, vjr = 0.0;
	double a_sum, b_sum, as, bs, ar, br;

	if (strcmp(scale, "old") == 0) {
		for (k = 1; k <= CatNum ; k++) {
			a_sum = k*D*a;
			b_sum = 0.0;
			for (l = 2; l <= k ; l++) {  /* b[1] = 0 */
				b_sum += b[l];
			}
			b_sum *= -D*a;

			/* new-to-old transformation */
			as = a_sum/S;
			bs = b_sum - (I/S)*a_sum;
			vjs += exp(as*theta + bs);  
		}
		a_sum = (CatID)*D*a;
		b_sum = 0.0;
		for (l = 2; l <= CatID ; l++) {
			b_sum += b[l];
		}
		b_sum *= -D*a;

		/* new-to-old transformation */
		as = a_sum/S;
		bs = b_sum - (I/S)*a_sum;
		return exp(as*theta + bs)/vjs;
	}
	else {
		for (k = 1; k <= CatNum ; k++) {
			a_sum = k*D*a;
			b_sum = 0.0;
			for (l = 2; l <= k ; l++) {
				b_sum += b[l];
			}
			b_sum *= -D*a;
			
			/* old-to-new transformation */
			ar = S*a_sum;
			br = b_sum + I*a_sum;
			vjr += exp(ar*theta + br);  
		}
		a_sum = (CatID)*D*a;
		b_sum = 0.0;
		for (l = 2; l <= CatID ; l++) {
			b_sum += b[l];
		}
		b_sum *= -D*a;

		/* old-to-new transformation */
		ar = S*a_sum;
		br = b_sum + I*a_sum;
		return exp(ar*theta + br)/vjr;  
	}
}

/***********************************************************************************/

double PdGPCOldOverS(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the GPC model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to S.
    
    Note:
	nb[2] for the actual first item-step parameter
	nb[CatNum] for the last item-step parameter
	is dealt with as a special case of the nominal response model

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double as, ps;
	double na_sum, as_ps_sum = 0.0;
	double ps_over_S;

	for (k = 1; k <= CatNum; k++) {
		na_sum = k*D*na;
		as = na_sum/S;      /* new-to-old transformation */
		ps = ProbGPC(CatNum, k, theta, D, na, nb, "old", S, I);
		as_ps_sum += as * ps;
	}
	   
	na_sum = CatID*D*na; /* for the category in question */
	as = na_sum/S;       /* new-to-old transformation */
	ps = ProbGPC(CatNum, CatID, theta, D, na, nb, "old", S, I);

	ps_over_S = -ps *((theta - I)/S)*(as - as_ps_sum);
	return ps_over_S;
}

/***********************************************************************************/

double PdGPCOldOverI(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the GPC model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to I.
    
    Note:
	nb[2] for the actual first item-step parameter
	nb[CatNum] for the last item-step parameter
	is dealt with as a special case of the nominal response model

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double as, ps;
	double na_sum, as_ps_sum = 0.0;
	double ps_over_I;

	for (k = 1; k <= CatNum; k++) {
		na_sum = k*D*na;
		as = na_sum/S;      /* new-to-old transformation */
		ps = ProbGPC(CatNum, k, theta, D, na, nb, "old", S, I);
		as_ps_sum += as * ps;
	}

	na_sum = CatID*D*na; /* for the category in question */
	as = na_sum/S;       /* new-to-old transformation */
	ps = ProbGPC(CatNum, CatID, theta, D, na, nb, "old", S, I);

	ps_over_I = -ps * (as - as_ps_sum);
	return ps_over_I;
}

/***********************************************************************************/

double PdGPCNewOverS(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the GPC model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to S.
    
    Note:
	ob[2] for the actual first item-step parameter
	ob[CatNum] for the last item-step parameter
	is dealt with as a special case of the nominal response model

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double pr;
	double oa_sum, oa_pr_sum = 0.0;
	double pr_over_S;

	for (k = 1; k <= CatNum; k++) {
		oa_sum = k*D*oa;
		pr = ProbGPC(CatNum, k, theta, D, oa, ob, "new", S, I);
		oa_pr_sum += oa_sum * pr;
	}
	oa_sum = CatID*D*oa; /* for the category in question */
	pr = ProbGPC(CatNum, CatID, theta, D, oa, ob, "new", S, I);
	pr_over_S = pr * theta *(oa_sum - oa_pr_sum);
	return pr_over_S;
}

/***********************************************************************************/

double PdGPCNewOverI(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the GPC model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to I.
    
    Note:
	ob[2] for the actual first item-step parameter
	ob[CatNum] for the last item-step parameter
	is dealt with as a special case of the nominal response model

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double pr;
	double oa_sum, oa_pr_sum = 0.0;
	double pr_over_I;

	for (k = 1; k <= CatNum ; k++) {
		oa_sum = k*D*oa;
		pr = ProbGPC(CatNum, k, theta, D, oa, ob, "new", S, I);
		oa_pr_sum += oa_sum * pr;
	}
	oa_sum = CatID*D*oa; /* for the category in question */
	pr = ProbGPC(CatNum, CatID, theta, D, oa, ob, "new", S, I);
	pr_over_I = pr * (oa_sum - oa_pr_sum);
	return pr_over_I;
}

/***********************************************************************************/

double ProbNRM(int CatNum, int CatID, double theta, double a[], double c[],
		char scale[], double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the nominal response model, calculate the category characteristic
    curve.
    
    Input
	CatNum: number of categories
	CatID: category response ID
	theta: ability value
	a[1..CatNum]: discrimination parameters
	c[1..CatNum]: intercept parameters
	scale: "old" or "new"
	S: slope of linear transformation
	I: intercept of linear transformation
	
	Note: No scaling constant
	
    Output
	Return either original or transformed probability.
	For the original probability, S = 1 and I = 0 with parameters on the
	reference scale.
	For the transfomred probability, S and I with parameters on the
	transformed scale.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double vjs = 0.0, vjr = 0.0;
	double as, cs, ar, cr;

	if(strcmp(scale, "old") == 0) {
		for (k = 1; k <= CatNum; k++) {
			as = a[k]/S;
			cs = c[k] - (I/S)*a[k];
			vjs += exp( as*theta + cs );
		}
		as = a[CatID]/S;
		cs = c[CatID] - (I/S)*a[CatID];
		return exp(as*theta + cs)/vjs;
	}
	else {
		for (k = 1; k <= CatNum; k++) {
			ar = S*a[k];
			cr = c[k] + I*a[k];
			vjr += exp( ar*theta + cr );
		}
		ar = S*a[CatID];
		cr = c[CatID] + I*a[CatID];
		return exp(ar*theta + cr)/vjr;
	}
}

/***********************************************************************************/

double PdNRMOldOverS(int CatNum, int CatID, double theta, double na[], double nc[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the NR model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to S.
    
    Note:
	na[1..CatNum] for the discrimination parameters
	nc[1..CatNum] for the intercept parameters
        No scaling constant

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double as, ps;
	double as_ps_sum = 0.0;
	double ps_over_S;

	for (k = 1; k <= CatNum; k++) {
		as = na[k]/S;
		ps = ProbNRM(CatNum, k, theta, na, nc, "old", S, I);
		as_ps_sum += as * ps;
	}   
	as = na[CatID]/S;
	ps = ProbNRM(CatNum, CatID, theta, na, nc, "old", S, I);
	ps_over_S = - ps*((theta - I)/S)*(as - as_ps_sum);
	return ps_over_S;
}

/***********************************************************************************/

double PdNRMOldOverI(int CatNum, int CatID, double theta, double na[], double nc[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the NR model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to I.
    
    Note:
	na[1..CatNum] for the discrimination parameters
	nc[1..CatNum] for the intercept parameters
        No scaling constant

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double as, ps;
	double as_ps_sum = 0.0;
	double ps_over_I;

	for (k = 1; k <= CatNum ; k++) {
		as = na[k]/S;
		ps = ProbNRM(CatNum, k, theta, na, nc, "old", S, I);
		as_ps_sum += as * ps;
	}   
		
	as = na[CatID]/S;
	ps = ProbNRM(CatNum, CatID, theta, na, nc, "old", S, I);
	ps_over_I = -ps * (as - as_ps_sum);
	return ps_over_I;
}

/***********************************************************************************/

double PdNRMNewOverS(int CatNum, int CatID, double theta, double oa[], double oc[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the NR model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to S.
    
    Note:
	oa[1..CatNum] for the discrimination parameters
	oc[1..CatNum] for the intercept parameters
        No scaling constant

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double pr;
	double oa_pr_sum = 0.0;
	double pr_over_S;

	for (k = 1; k <= CatNum ; k++) {
		pr = ProbNRM(CatNum, k, theta, oa, oc, "new", S, I);
		oa_pr_sum += oa[k] * pr;
	}   
	pr = ProbNRM(CatNum, CatID, theta, oa, oc, "new", S, I);
	pr_over_S = pr*theta*(oa[CatID] - oa_pr_sum);
	return pr_over_S;
}

/***********************************************************************************/

double PdNRMNewOverI(int CatNum, int CatID, double theta, double oa[], double oc[],
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the NR model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to I.
    
    Note:
	oa[1..CatNum] for the discrimination parameters
	oc[1..CatNum] for the intercept parameters
        No scaling constant

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int k;
	double pr;
	double oa_pr_sum = 0.0;
	double pr_over_I;

	for (k = 1; k <= CatNum ; k++) {
		pr = ProbNRM(CatNum, k, theta, oa, oc, "new", S, I);
		oa_pr_sum += oa[k] * pr;
	}   
		
	pr = ProbNRM(CatNum, CatID, theta, oa, oc, "new", S, I);
	pr_over_I = pr*(oa[CatID] - oa_pr_sum);
	return pr_over_I;
}          

/***********************************************************************************/

void StHaebara(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		enum symmetry SYM, enum OnOff FuncStd, double S0, double I0,
		double *slope, double *intercept)
/*------------------------------------------------------------------------------
  Functionality:
    Estimate the slope (A) and intercept (B) of the new-to-old transformation
    by using the Haebara method.
    
    Input:
	Handle: A pointer to control solutions
	ComItem: A pointer to an array of the CommonItemSpec structure
	SYM: Symmetric option (old_scale, new_scale, or symmetric)
	FuncStd: Function standardization option (on or off)
	S0: A starting value for the slope
	I0: A starting value of the intercept
    Output:    
	*slope: estimate of A
	*intercept: estimate of B

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int n, iter;
	double x[3];
	double ftol=0.0000000001, fret;

	ContHandle = Handle;
	ContComItem = ComItem;
	ContSym = SYM;
	ContFuncStd = FuncStd;

	x[1] = S0;
	x[2] = I0;
	n = 2;
	/* ===== Before version 1.0 update =====
	dfpmin(x, n, ftol, &iter, &fret, FuncHaebara, GradHaebara);
	*/
	er_dfpmin(x, n, ftol, &iter, &fret, FuncHaebara, GradHaebara);
	/* ===== End of version 1.0 update ===== */ 
	*slope = x[1];
	*intercept = x[2];
	return;
}

/***********************************************************************************/

double FuncHaebara(double x[])
/*------------------------------------------------------------------------------
  Functionality:
    calculate the value of the Haebara criterion function.
    
    Input:
	x: the point of slope and intercept
    Output:    
	the value of the criterion function

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int i, j, k;
	int cat_sum = 0;
	double w1_sum = 0.0, w2_sum = 0.0;
	double sym_f1, sym_f2;
	double theta, th_weight;
	double p_origi, p_trans;
	double func1, func2;
	double func1_sum = 0.0, func2_sum = 0.0;
 
	/* Q1: old scale */
	for (i = 1; i <= ContHandle->OldThetaNum; i++) {
		theta = ContHandle->OldThetaValues[i];
		th_weight = ContHandle->OldThetaWeights[i];
		w1_sum += th_weight;

		func1 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				p_origi = ProbOld(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbOld(&ContComItem[j], k, theta, off, x[1], x[2]);
				func1 += (p_origi - p_trans)*(p_origi - p_trans);
			}
		}
		func1_sum += func1 * th_weight;
	}

	/* Q2: new scale */
	for (i = 1; i <= ContHandle->NewThetaNum; i++) {
		theta = ContHandle->NewThetaValues[i];
		th_weight = ContHandle->NewThetaWeights[i];
		w2_sum += th_weight;
		
		func2 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				p_origi = ProbNew(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbNew(&ContComItem[j], k, theta, off, x[1], x[2]);
				func2 += (p_origi - p_trans)*(p_origi - p_trans);
			}
			if (i == 1) cat_sum += ContComItem[j].CatNum;
		}
		func2_sum += func2 * th_weight;
	}

	/* symmetric or non-symmetric setting*/
	if (ContSym == symmetric) {sym_f1 = 1.0; sym_f2 = 1.0;}
	else {
		if (ContSym == old_scale) {sym_f1 = 1.0; sym_f2 = 0.0;}
		else                      {sym_f1 = 0.0; sym_f2 = 1.0;}
	}

	/* function standardization */
	if (ContFuncStd == on) {
		return (sym_f1*func1_sum/(cat_sum*w1_sum) + sym_f2*func2_sum/(cat_sum*w2_sum));
	}
	else {
		return (sym_f1*func1_sum + sym_f2*func2_sum);
	}
}

/***********************************************************************************/

void GradHaebara(double x[], double grad[])
/*------------------------------------------------------------------------------
  Functionality:
    calculate the partial derivatives of the Haebara criterion function with
    respect to slope (S or A) and intercept (I or B).
    
    Input:
	x: the point of slope and intercept
    Output:    
	the gradient of the criterion function w.r.t. A and B.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int i, j, k;
	int cat_sum = 0;
	double w1_sum = 0.0, w2_sum = 0.0;
	double sym_f1, sym_f2;
	double theta, th_weight;
	double p_origi, p_trans, ps, pi;
	double ps_f1, pi_f1, ps_f2, pi_f2;
	double ps_f1_sum = 0.0, pi_f1_sum = 0.0;
	double ps_f2_sum = 0.0, pi_f2_sum = 0.0;

	/* Q1: old scale */
	for (i = 1; i <= ContHandle->OldThetaNum; i++) {
		theta = ContHandle->OldThetaValues[i];
		th_weight = ContHandle->OldThetaWeights[i];
		w1_sum += th_weight;

		ps_f1 = 0.0; pi_f1 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				p_origi = ProbOld(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbOld(&ContComItem[j], k, theta, off, x[1], x[2]);
				ps = PdOldOverS(&ContComItem[j], k, theta, x[1], x[2]);
				pi = PdOldOverI(&ContComItem[j], k, theta, x[1], x[2]);

				ps_f1 += (p_origi - p_trans)*ps;
				pi_f1 += (p_origi - p_trans)*pi;
			}
			if (i == 1) cat_sum += ContComItem[j].CatNum;
		}
		ps_f1_sum += ps_f1 * th_weight;
		pi_f1_sum += pi_f1 * th_weight;
	}

	/* Q2: new scale */

	for (i = 1; i <= ContHandle->NewThetaNum; i++) {
		theta = ContHandle->NewThetaValues[i];
		th_weight = ContHandle->NewThetaWeights[i];
		w2_sum += th_weight;

		ps_f2 = 0.0; pi_f2 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				p_origi = ProbNew(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbNew(&ContComItem[j], k, theta, off, x[1], x[2]);
				ps = PdNewOverS(&ContComItem[j], k, theta, x[1], x[2]);
				pi = PdNewOverI(&ContComItem[j], k, theta, x[1], x[2]);

				ps_f2 += (p_origi - p_trans)*ps;
				pi_f2 += (p_origi - p_trans)*pi;
			}
		}
		ps_f2_sum += ps_f2 * th_weight;
		pi_f2_sum += pi_f2 * th_weight;
	}

	/* symmetric or non-symmetric setting*/
	if (ContSym == symmetric) {sym_f1 = 1.0; sym_f2 = 1.0;}
	else {
		if (ContSym == old_scale) {sym_f1 = 1.0; sym_f2 = 0.0;}
		else                      {sym_f1 = 0.0; sym_f2 = 1.0;}
	}

	/* function standardization */
	if (ContFuncStd == on) {
		grad[1] = -2.0*(sym_f1*ps_f1_sum/(cat_sum*w1_sum) + sym_f2*ps_f2_sum/(cat_sum*w2_sum));
		grad[2] = -2.0*(sym_f1*pi_f1_sum/(cat_sum*w1_sum) + sym_f2*pi_f2_sum/(cat_sum*w2_sum));
	}
	else {
		grad[1] = -2.0*(sym_f1*ps_f1_sum + sym_f2*ps_f2_sum);
		grad[2] = -2.0*(sym_f1*pi_f1_sum + sym_f2*pi_f2_sum);

	}
	return;
}

/***********************************************************************************/

void StMeanMean(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		double *slope, double *intercept)
/*------------------------------------------------------------------------------
  Functionality:
    Estimate the slope (A) and intercept (B) of the new-to-old transformation
    by using the Mean/Mean method.
    
    Input:
	Handle: A pointer to an array of the CommonItemSpec structure
	ComItem: A pointer to the CommonItemSpec structure
    Output:    
	*slope: estimate of A
	*intercept: estimate of B

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int j, k, a_num=0, b_num=0, index_a=0, index_b=0;
	double new_mu_a=0.0, old_mu_a=0.0;
	double new_mu_b=0.0, old_mu_b=0.0;
	double *na_vec, *oa_vec, *nb_vec, *ob_vec;

	/* counting valid a- and b-parameters by model */

	for(j=1; j <= Handle->ComItemNum; j++) {
		if(ComItem[j].model==l3 || ComItem[j].model==gr ||ComItem[j].model==pc) {
			a_num++;
			b_num += (ComItem[j].CatNum-1);
		}
		else {
			a_num += ComItem[j].CatNum;
			b_num += ComItem[j].CatNum;
		}
	}

	/* memory allocation */
	na_vec = (double *) malloc((a_num+1)*sizeof(double));
	if (!na_vec) 
		runerror("memory allocation failure 1 in StMeanMean()\n");

	oa_vec = (double *) malloc((a_num+1)*sizeof(double));
	if (!oa_vec) 
		runerror("memory allocation failure 2 in StMeanMean()\n");

	nb_vec = (double *) malloc((b_num+1)*sizeof(double));
	if (!nb_vec)
		runerror("memory allocation failure 3 in StMeanMean()\n");

	ob_vec = (double *) malloc((b_num+1)*sizeof(double));
	if (!ob_vec) 
		runerror("memory allocation failure 4 in StMeanMean()\n");


	/* copying valid a- and b-parameters by model into nb_vec and ob_vec */
	for(j=1; j <= Handle->ComItemNum; j++) {
		if(ComItem[j].model==l3 || ComItem[j].model==gr ||ComItem[j].model==pc) {
			index_a++;
			na_vec[index_a] = ComItem[j].Na[2];
			oa_vec[index_a] = ComItem[j].Oa[2];
			for(k=2; k <= ComItem[j].CatNum; k++) {
				index_b++;
				nb_vec[index_b] = ComItem[j].Nb[k];
				ob_vec[index_b] = ComItem[j].Ob[k];
			}
		}
	 	else {
			for(k=1; k <= ComItem[j].CatNum; k++) {
				index_a++;
				index_b++;
				na_vec[index_a] = ComItem[j].Na[k];
				oa_vec[index_a] = ComItem[j].Oa[k];
				nb_vec[index_b] = -(ComItem[j].Nc[k]/ComItem[j].Na[k]);
				ob_vec[index_b] = -(ComItem[j].Oc[k]/ComItem[j].Oa[k]);
			}
		}
	}

	for (j=1; j <= a_num; j++) {
		new_mu_a += na_vec[j];
		old_mu_a += oa_vec[j];
	}
	for (j=1; j <= b_num; j++) {
		new_mu_b += nb_vec[j];
		old_mu_b += ob_vec[j];
	}
	if (a_num != 0 && b_num != 0) {
		new_mu_a /= a_num;
		old_mu_a /= a_num;
		new_mu_b /= b_num;
		old_mu_b /= b_num;

		*slope = new_mu_a/old_mu_a;
		*intercept = old_mu_b - (*slope)*new_mu_b;
	}
	else {
		*slope = 1.0;
		*intercept = 0.0;
	}

	free( (void *) na_vec);
	free( (void *) oa_vec);
	free( (void *) nb_vec);
	free( (void *) ob_vec);
	return;
}

/***********************************************************************************/

void StMeanSigma(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		double *slope, double *intercept)
/*------------------------------------------------------------------------------
  Functionality:
    Estimate the slope (A) and intercept (B) of the new-to-old transformation
    by using the Mean/Sigma method.
    
    Input:
	Handle: A pointer to an array of the CommonItemSpec structure
	ComItem: A pointer to the CommonItemSpec structure
    Output:    
	*slope: estimate of A
	*intercept: estimate of B

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int j, k, b_num=0, index=0;
	double new_mu_b=0.0, old_mu_b=0.0;
	double new_si_b=0.0, old_si_b=0.0;
	double *nb_vec, *ob_vec;

	/* counting valid b-parameters by model */

	for(j=1; j <= Handle->ComItemNum; j++) {
		if(ComItem[j].model==l3 || ComItem[j].model==pc ||ComItem[j].model==gr) {
			b_num += (ComItem[j].CatNum-1);
		}
		else b_num += ComItem[j].CatNum;
	}

	/* memory allocation */
	nb_vec = (double *) malloc((b_num+1)*sizeof(double));
	if (!nb_vec) 
		runerror("memory allocation failure 1 in StMeanSigma()\n");

	ob_vec = (double *) malloc((b_num+1)*sizeof(double));
	if (!ob_vec) 
		runerror("memory allocation failure 2 in StMeanSigma()\n");


	/* copying valid b-parameters by model into nb_vec and ob_vec */
	for(j=1; j <= Handle->ComItemNum; j++) {
		if(ComItem[j].model==l3 || ComItem[j].model==pc ||ComItem[j].model==gr) {
			for(k=2; k <= ComItem[j].CatNum; k++) {
				index++;
				nb_vec[index] = ComItem[j].Nb[k];
				ob_vec[index] = ComItem[j].Ob[k];
			}
		}
		else {
			for(k=1; k <= ComItem[j].CatNum; k++) {
				index++;
				nb_vec[index] = -(ComItem[j].Nc[k]/ComItem[j].Na[k]);
				ob_vec[index] = -(ComItem[j].Oc[k]/ComItem[j].Oa[k]);
			}
		}
	}

	for (j=1; j <= b_num; j++) {
		new_mu_b += nb_vec[j];
		old_mu_b += ob_vec[j];
		new_si_b += nb_vec[j]*nb_vec[j];
		old_si_b += ob_vec[j]*ob_vec[j];
	}
	if (b_num != 0) {
		new_mu_b /= b_num;
		old_mu_b /= b_num;
		new_si_b = sqrt(new_si_b/b_num - new_mu_b*new_mu_b);
		old_si_b = sqrt(old_si_b/b_num - old_mu_b*old_mu_b);

		*slope = old_si_b/new_si_b;
		*intercept = old_mu_b - (*slope)*new_mu_b;
	}
	else {
		*slope = 1.0;
		*intercept = 0.0;
	}

	free( (void *) nb_vec);
	free( (void *) ob_vec);
	return;
}

/***********************************************************************************/

void StStockingLord(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		enum symmetry SYM, enum OnOff FuncStd, double S0, double I0,
		double *slope, double *intercept)
/*------------------------------------------------------------------------------
  Functionality:
    Estimate the slope (A) and intercept (B) of the new-to-old transformation
    by using the Stocking-Lord method.
    
    Input:
	Handle: A pointer to control solutions
	ComItem: A pointer to an array of the CommonItemSpec structure
	SYM: Symmetric option (old_scale, new_scale, symmetric)
	FuncStd: Function standardization option (on or off)
	S0: A starting value for the slope
	I0: A starting value of the intercept
    Output:    
	*slope: estimate of A
	*intercept: estimate of B

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int n, iter;
	double x[3];
	double ftol=0.0000000001, fret;

	ContHandle = Handle;
	ContComItem = ComItem;
	ContSym = SYM;
	ContFuncStd = FuncStd;

	x[1] = S0;
	x[2] = I0;
	n = 2;

	/* ===== Before version 1.0 update =====
	dfpmin(x, n, ftol, &iter, &fret, FuncStockingLord, GradStockingLord);
	*/
	er_dfpmin(x, n, ftol, &iter, &fret, FuncStockingLord, GradStockingLord); 
	/* ===== End of version 1.0 update ===== */ 

	*slope = x[1];
	*intercept = x[2];
	return;
}

/***********************************************************************************/

double FuncStockingLord(double x[])
/*------------------------------------------------------------------------------
  Functionality:
    calculate the value of the Stocking-Lord criterion function.
    
    Input:
	x: the point of slope and intercept
    Output:    
	the value of the criterion function

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int i, j, k;
	double w1_sum = 0.0, w2_sum = 0.0;
	double sym_f1, sym_f2;
	double theta, th_weight;
	double p_origi, p_trans, Ujk;
	double func1, func2;
	double func1_sum = 0.0, func2_sum = 0.0;

	/* F1: old scale */
	for (i = 1; i <= ContHandle->OldThetaNum; i++) {
		theta = ContHandle->OldThetaValues[i];
		th_weight = ContHandle->OldThetaWeights[i];
		w1_sum += th_weight;

		func1 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				Ujk = ContComItem[j].ScoreFunc[k];
				p_origi = ProbOld(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbOld(&ContComItem[j], k, theta, off, x[1], x[2]);

				func1 += Ujk*(p_origi - p_trans);
			}
		}
		func1_sum += (func1 * func1) * th_weight;
	}

	/* F2: new scale */
	for (i = 1; i <= ContHandle->NewThetaNum; i++) {
		theta = ContHandle->NewThetaValues[i];
		th_weight = ContHandle->NewThetaWeights[i];
		w2_sum += th_weight;

		func2 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				Ujk = ContComItem[j].ScoreFunc[k];
				p_origi = ProbNew(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbNew(&ContComItem[j], k, theta, off, x[1], x[2]);

				func2 += Ujk*(p_origi - p_trans);
			}
		}
		func2_sum += (func2 * func2) * th_weight;
	}

	/* symmetric or non-symmetric setting*/
	if (ContSym == symmetric) {sym_f1 = 1.0; sym_f2 = 1.0;}
	else {
		if (ContSym == old_scale) {sym_f1 = 1.0; sym_f2 = 0.0;}
		else                      {sym_f1 = 0.0; sym_f2 = 1.0;}
	}

	/* function standardization */
	if (ContFuncStd == on) {
		return (sym_f1*func1_sum/w1_sum + sym_f2*func2_sum/w2_sum);
	}
	else {
		return (sym_f1*func1_sum + sym_f2*func2_sum);
	}
}

/***********************************************************************************/

void GradStockingLord(double x[], double grad[])
/*------------------------------------------------------------------------------
  Functionality:
    calculate the partial derivatives of the Stocking-Lord criterion function
    with respect to slope (S or A) and intercept (I or B).
    
    Input:
	x: the point of slope and intercept
    Output:    
	the gradient of the criterion function w.r.t. A and B.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int i, j, k;
	double w1_sum = 0.0, w2_sum = 0.0;
	double sym_f1, sym_f2;

	double theta, th_weight;
	double p_origi, p_trans, ps, pi, Ujk;
	double func1, func2, ps_sum, pi_sum;
	double ps_f1_sum = 0.0, pi_f1_sum = 0.0;
	double ps_f2_sum = 0.0, pi_f2_sum = 0.0;


	/* F1: old scale */
	for (i = 1; i <= ContHandle->OldThetaNum; i++) {
		theta = ContHandle->OldThetaValues[i];
		th_weight = ContHandle->OldThetaWeights[i];
		w1_sum += th_weight;

		func1 = 0.0;
		ps_sum = 0.0; pi_sum = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				Ujk = ContComItem[j].ScoreFunc[k];
				p_origi = ProbOld(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbOld(&ContComItem[j], k, theta, off, x[1], x[2]);
				ps = PdOldOverS(&ContComItem[j], k, theta, x[1], x[2]);
				pi = PdOldOverI(&ContComItem[j], k, theta, x[1], x[2]);

				func1 += Ujk * (p_origi - p_trans);
				ps_sum += Ujk*ps;
				pi_sum += Ujk*pi;
			}
		}
		ps_f1_sum += func1 * ps_sum * th_weight;
		pi_f1_sum += func1 * pi_sum * th_weight;
	}

	/* F2: new scale */
	for (i = 1; i <= ContHandle->NewThetaNum; i++) {
		theta = ContHandle->NewThetaValues[i];
		th_weight = ContHandle->NewThetaWeights[i];
		w2_sum += th_weight;

		func2 = 0.0;
		ps_sum = 0.0; pi_sum = 0.0;

		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				Ujk = ContComItem[j].ScoreFunc[k];
				p_origi = ProbNew(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbNew(&ContComItem[j], k, theta, off, x[1], x[2]);
				ps = PdNewOverS(&ContComItem[j], k, theta, x[1], x[2]);
				pi = PdNewOverI(&ContComItem[j], k, theta, x[1], x[2]);

				func2 += Ujk * (p_origi - p_trans);
				ps_sum += Ujk*ps;
				pi_sum += Ujk*pi;
			}
		}
		ps_f2_sum += func2 * ps_sum * th_weight;
		pi_f2_sum += func2 * pi_sum * th_weight;
	}

	/* symmetric or non-symmetric setting*/
	if (ContSym == symmetric) {sym_f1 = 1.0; sym_f2 = 1.0;}
	else {
		if (ContSym == old_scale) {sym_f1 = 1.0; sym_f2 = 0.0;}
		else                      {sym_f1 = 0.0; sym_f2 = 1.0;}
	}

	/* function standardization */
	if (ContFuncStd == on) {
		grad[1] = -2.0*(sym_f1*ps_f1_sum/w1_sum + sym_f2*ps_f2_sum/w2_sum);
		grad[2] = -2.0*(sym_f1*pi_f1_sum/w1_sum + sym_f2*pi_f2_sum/w2_sum);
	}
	else {
		grad[1] = -2.0*(sym_f1*ps_f1_sum + sym_f2*ps_f2_sum);
		grad[2] = -2.0*(sym_f1*pi_f1_sum + sym_f2*pi_f2_sum);
	}
	return;
}

/***********************************************************************************/

void ThetaInfoRead(FILE *inf, const char *oldOrnew, struct IRTstControl *Handle)
/*------------------------------------------------------------------------------
  Functionality:
    Read pieces of information for ability distributions of the new and old
    groups from an input file, depending on the oldOrnew string, which must be
    "new" or "old."
    Assume that the input file is opened but is not read in yet at all.
    While reading the pairs information, the member pointers of IRTstControl are
    set up by dynamic memory allocation. In addition, some members of
    an IRTstControl structure referenced by the pointer Handle are initialized
    to their respective values that have been read.
    
    Input:
      inf       A file pointer; the file must be opened for input.                       
      oldOrnew  A string indicating that the ability distribution is for the old
                or new group. The string must be either "new" or "old"
      Handle    A pointer to a variable of the IRTstControl structure, which is
                defined in the header file IRTst.h 

      The input file must be prepared according to the following format:
      Line 1    Number of the ability points for the new or old group.
      Line 2+   Pairs of ability points and weights.
                Each pair consists of an ability point and its weight.

                Suppose, for example, that there are five ability points to be used,
                either on the new or old scale. Pairs of points and weights for 
                ability, then, must be typed as follows, with the first line being
                for the number of the pairs:

                5
                -2.0 0.05
                -1.0 0.25
                 0.0 0.40
                 1.0 0.25
                 2.0 0.05
              
    While reading the pairs information, the member pointers of IRTstControl are
    set by dynamic memory allocation.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int i, ThetaNum;

	fscanf(inf, "%d",   &ThetaNum); /* reading the number of ability points */
	if (strcmp(oldOrnew, "old") == 0) Handle->OldThetaNum = ThetaNum;
	else                              Handle->NewThetaNum = ThetaNum;

	/* memory allocation */
	if (strcmp(oldOrnew, "old") == 0) {
	        Handle->OldThetaValues = (double *) malloc((ThetaNum+1)*sizeof(double));
		if (!Handle->OldThetaValues) 
			runerror("memory allocation failure 1 in ThetaInfoRead()\n");

	        Handle->OldThetaWeights = (double *) malloc((ThetaNum+1)*sizeof(double));
		if (!Handle->OldThetaWeights) 
			runerror("memory allocation failure 2 in ThetaInfoRead()\n");

	}
	else {
	        Handle->NewThetaValues = (double *) malloc((ThetaNum+1)*sizeof(double));
		if (!Handle->NewThetaValues) 
			runerror("memory allocation failure 3 in ThetaInfoRead()\n");

	        Handle->NewThetaWeights = (double *) malloc((ThetaNum+1)*sizeof(double));
		if (!Handle->NewThetaWeights) 
			runerror("memory allocation failure 4 in ThetaInfoRead()\n");

	}
        
	for (i=1; i <= ThetaNum; i++) {
		if (strcmp(oldOrnew, "old") == 0) {
	        	fscanf(inf, "%lf", &(Handle->OldThetaValues[i]));
	        	fscanf(inf, "%lf", &(Handle->OldThetaWeights[i]));
		}
		else {
	        	fscanf(inf, "%lf", &(Handle->NewThetaValues[i]));
	        	fscanf(inf, "%lf", &(Handle->NewThetaWeights[i]));
		}
	}
	return;
}

/***********************************************************************************/

void StItemDeAlloc(struct ItemSpec *Item, const char *oldOrnew, struct IRTstControl *Handle)
/*------------------------------------------------------------------------------
  Functionality:
    Deallocate memory given to the pointer (Item) to an ItemSpec structure,
    depending on the scale ("new" or "old").
    The memory given to its members is also deallocated.
    
    Input:
      Item      A pointer to an array of the ItemSpec structure
      oldOrnew  A string indicating that the array is about either the old or
                new form
      Handle    A pointer to a variable of the IRTstControl structure

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int j, ItemNum;
	if (strcmp(oldOrnew, "old")==0) ItemNum = Handle->OldItemNum;
	else                            ItemNum = Handle->NewItemNum;
	
	for (j=1; j <= ItemNum; j++) {
		free( (void *) Item[j].ScoreFunc);
		free( (void *) Item[j].a);
		free( (void *) Item[j].b);
		free( (void *) Item[j].c);
		free( (void *) Item[j].d);
	}

	free( (void *) Item);
	
	return;
}

/***********************************************************************************/

void StComItemDeAlloc(struct CommonItemSpec *ComItem, struct IRTstControl *Handle)
/*------------------------------------------------------------------------------
  Functionality:
    Deallocate memory given to the pointer (ComItem) to a CommonItemSpec
    structure.
    The memory given to its members is also deallocated.
    
    Input:
      ComItem   A pointer to an array of the CommonItemSpec structure
      Handle    A pointer to a variable of the IRTstControl structure

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	int j;
	for (j=1; j <= Handle->ComItemNum; j++) {
		free( (void *) ComItem[j].ScoreFunc);
		free( (void *) ComItem[j].Na);
		free( (void *) ComItem[j].Nb);
		free( (void *) ComItem[j].Nc);
		free( (void *) ComItem[j].Nd);
		free( (void *) ComItem[j].Oa);
		free( (void *) ComItem[j].Ob);
		free( (void *) ComItem[j].Oc);
		free( (void *) ComItem[j].Od);
	}
	free( (void *) ComItem);
	
	return;
}

/***********************************************************************************/

void StContDeAlloc(struct IRTstControl *Handle)
/*------------------------------------------------------------------------------
  Functionality:
    Deallocate memory given to the members of the pointer Handle, which points
    to an IRTstControl structure.
    
    Input:
      Handle    A pointer to a variable of the IRTstControl structure

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	free( (void *) Handle->NewThetaValues);
	free( (void *) Handle->NewThetaWeights);
	free( (void *) Handle->OldThetaValues);
	free( (void *) Handle->OldThetaWeights);

	return;
}

/***********************************************************************************/

 void Wrapper_IRTst(FILE *outf, char tt[],
     char ItemNewFile[], char ItemOldFile[], char ItemCommonFile[], 
     char DistNewFile[], char DistOldFile[], 
     int HA, enum symmetry HAsym, enum OnOff HAfs, double HAs, double HAi,
     int SL, enum symmetry SLsym, enum OnOff SLfs, double SLs, double SLi,
     char ST[], int PrintFiles)         
/*
  Wrapper for IRT scale-score transformation.  As distinct from other Wrapper
  functions in Equating Recipes, this Wrapper function does its own printing
  
  Input:
  
    outf:             pointer to output file
    tt                user-supplied title
    ItemNewFile[]:    name of file containing item parameters for new form X 
    ItemOldFile[]:    name of file containing item parameters for old form Y
    ItemCommonFile[]: name of file containing id's for common items
    DistNewFile[]:    name of file containing quadrature ability scale for new form X
    DistOldFile[]:    name of file containing quadrature ability scale for old form Y 
    HA:               Haebara results: 1 --> compute results; 0 --> no results
	HAsym             Haebara symmetric option (old_scale, new_scale, symmetric)
	HAfs:             Haebara function standardization option (on or off)
	HAs:              Haebara starting value for the slope
	HAi:              Haebara starting value of the intercept
    SL:               Stocking-Lord results: 1 --> compute results; 0 --> no results
	SLsym             Stocking-Lord symmetric option (old_scale, new_scale, symmetric)
	SLfs:             Stocking-Lord function standardization option (on or off)
	SLs:              Stocking-Lord starting value for the slope
	SLi:              Stocking-Lord starting value for the intercept 
    ST:               method to be used for scale transformation
                        ("MS" --> mean/sigma; "MM" --> mean/mean;
                         "HA" --> Haebara: "SL" --> Stocking-Lord; 
                         "NL" --> no scale transformation)
    PrintFiles: 0 --> don't print input item parameter and quadrature files
                1 --> print input item parameter and quadrature files

  Note: If either DistNewFile[] or DistOldFile[] is set to NULL,
        then only mean/mean and mean/sigma are computed
                            
  NOTE:  No capability built in for bootstrapping

  NOTE:  Structures are NOT saved after function returns
                                            
  Function calls other than C or NR utilities:                    
    ItemInfoRead()
    ItemPairsRead()
    ThetaInfoRead()
    StMeanSigma()
    StMeanMean()
    StHaebara()
	StStockingLord()
    ScaleTransform()
    StComItemDeAlloc()
    StItemDeAlloc()
    StItemDeAlloc()
    StContDeAlloc()
                                               
  Author: R. L. Brennan 
  Date of last revision 9/25/08

*/
{ 
  FILE *inf;
  struct ItemSpec *NewItems, *OldItems;
  struct CommonItemSpec *ComItems;
  struct IRTstControl StInfo; 
  double MSa, MSb, MMa, MMb, HAa, HAb, SLa, SLb, a, b; 
  char ItemNewFileST[103],    /* file name for transformed new item parameters */
       DistNewFileST[103]; /* file name for transformed new group ability dist */
  char line1[1000],                            /* first line; \n replace by \t */
       line2[1000];                                             /* second line */                           

  /* read the input for the new form items */

  inf = fopen(ItemNewFile, "r"); 
  if(!inf) runerror("can't open file for new form items \n");
  NewItems = ItemInfoRead(inf, "new", &StInfo);
  fclose(inf);

  /* read the input for the old form items */

  inf = fopen(ItemOldFile, "r");
  if(!inf) runerror("can't open file for old form items \n");
  OldItems = ItemInfoRead(inf, "old", &StInfo);
  fclose(inf);  

  /* read the input for the common items */

  inf = fopen(ItemCommonFile, "r");
  if(!inf) runerror("can't open file for common items\n");
  ComItems = ItemPairsRead(inf, NewItems, OldItems, &StInfo);
  fclose(inf);

  if(strcmp(DistNewFile,"NULL") != 0 && strcmp(DistOldFile,"NULL") != 0){

    /* read the input for the new group's ability distribution */

    inf = fopen(DistNewFile, "r");
    if(!inf) runerror("can't open file for new form ability distribution \n");
    ThetaInfoRead(inf, "new", &StInfo);
    fclose(inf);

    /* read the input for the old group's ability distribution */

    inf = fopen(DistOldFile, "r");
    if(!inf) runerror("can't open file for old form ability distribution \n");
    ThetaInfoRead(inf, "old", &StInfo);
    fclose(inf);
  }

  /* initial print statements */

  fprintf(outf,"\n\n%s\n\n",tt);
  fprintf(outf,"IRT Scale Transformation\n\n");

  fprintf(outf,"Input files:\n\n");

  fprintf(outf,"   %s (item parameters for new form X)\n",ItemNewFile);
  if(PrintFiles == 1) Print_file(ItemNewFile,outf);

  fprintf(outf,"   %s (item parameters for new form Y)\n",ItemOldFile);
  if(PrintFiles == 1) Print_file(ItemOldFile,outf);

  fprintf(outf,"   %s (pairs of ID's for common items)\n",ItemCommonFile);
  if(PrintFiles == 1) Print_file(ItemCommonFile,outf);

  if(strcmp(DistNewFile,"NULL") != 0 && strcmp(DistOldFile,"NULL") != 0){
    fprintf(outf,"   %s (ability distribution for new form X)\n",DistNewFile);
    if(PrintFiles == 1) Print_file(DistNewFile,outf);
    fprintf(outf,"   %s (ability distribution for old form Y)\n\n",DistOldFile);
    if(PrintFiles == 1) Print_file(DistOldFile,outf);
  }

  if(HA == 1){
    fprintf(outf,"Haebara symmetric option = %s\n",
             (HAsym==0 ? "old_scale" :(HAsym==1 ? "new_scale" : "symmetric")));
    fprintf(outf,"Haebara function standardization option = %s\n",
             HAfs ? "on" : "off"); 
    fprintf(outf,"Haebara starting value for slope (A) = %9.5f\n",HAs);
    fprintf(outf,"Haebara starting value for intercept (B) = %9.5f\n\n",HAi);
  }

  if(SL == 1){
    fprintf(outf,"Stocking-Lord symmetric option = %s\n",
             (SLsym==0 ? "old_scale" :(SLsym==1 ? "new_scale" : "symmetric")));
    fprintf(outf,"Stocking-Lord function standardization option = %s\n",
             SLfs ? "on" : "off"); 
    fprintf(outf,"Stocking-Lord starting value for slope (A) = %9.5f\n",SLs);
    fprintf(outf,"Stocking-Lord starting value for intercept (B) = %9.5f\n\n",SLi);
  }

  fprintf(outf,"------------------------------------------------\n");
  fprintf(outf,"                         A=slope     B=intercept\n");
  
  /* mean/sigma method */

  StMeanSigma(&StInfo, ComItems, &MSa, &MSb);
  fprintf(outf,"Mean/Sigma Method:    %10.5f      %10.5f\n", MSa, MSb);
	
  /* mean/mean method */

  StMeanMean(&StInfo, ComItems, &MMa, &MMb);
  fprintf(outf,"Mean/Mean Method:     %10.5f      %10.5f\n", MMa, MMb);

  if(HA != 1 && SL != 1)
    fprintf(outf,"------------------------------------------------\n");

  /* Haebara method */

  if(HA == 1){                                 /* Haebara method requested */
    if(strcmp(DistNewFile,"NULL") == 0 || strcmp(DistOldFile,"NULL") == 0)
       runerror("Haebara results cannot be computed because at least\n"
                "one of the ability distribution files is not present");
    StHaebara(&StInfo, ComItems, HAsym, HAfs, HAs, HAi, &HAa, &HAb);
    fprintf(outf,"Haebara Method:       %10.5f      %10.5f\n", HAa, HAb);
  }

  if(HA == 1 && SL != 1)
    fprintf(outf,"------------------------------------------------\n");

  /* Stocking-Lord method */

  if(SL == 1){                            /* Stocking-Lord method requested */
    if(strcmp(DistNewFile,"NULL") == 0 || strcmp(DistOldFile,"NULL") == 0)
      runerror("Stocking-Lord results cannot be computed because at least\n"
               "one of the ability distribution files is not present");
	StStockingLord(&StInfo, ComItems, SLsym, SLfs, SLs, SLi, &SLa, &SLb);
	fprintf(outf,"Stocking-Lord Method: %10.5f      %10.5f\n", SLa, SLb);
  }
  if(HA == 1 && SL == 1)
    fprintf(outf,"------------------------------------------------\n");

  /* scale transformation */ 
   
  if(strcmp(ST,"NL") != 0){

    strcat(strcpy(ItemNewFileST, ItemNewFile),".ST");
    strcat(strcpy(DistNewFileST, DistNewFile),".ST");

    if(strcmp(ST,"MS") == 0)      {a = MSa; b = MSb;}
    else if(strcmp(ST,"MM") == 0) {a = MMa; b = MMb;}
    else if(strcmp(ST,"HA") == 0) {a = HAa; b = HAb;}
    else if(strcmp(ST,"SL") == 0) {a = SLa; b = SLb;}
    else runerror("invalid scale transformation method"); 

    ScaleTransform(ItemNewFileST, DistNewFileST, a, b, NewItems, &StInfo);

    /* print Transformed Item Paramter Estimates for New Form */

    fprintf(outf,"\nTransformed Item Paramter Estimates for New Form\n"
                 "Based on %s Method\n"
                 "Results stored in file named: %s\n\n",
                 (!strcmp(ST,"MS") ? "Mean/Sigma" :
                 (!strcmp(ST,"MM") ? "Mean/Mean" : 
                 (!strcmp(ST,"HA") ? "Haebara" : "Stocking-Lord"))),  
           ItemNewFileST);

    inf = fopen(ItemNewFileST,"r");
    fgets(line1,998,inf);                                 /* skip the first line */
    while(fgets(line1,998,inf) != NULL){
      line1[strlen(line1)-1] = '\t';                      /*  replace \n with \t */
      fprintf(outf,"%s",strcat(line1,fgets(line2,998,inf)));
    }  
    fclose(inf);
    fprintf(outf,"------------------------------------------------\n");

    /* print Transformed Ability Distribution for New Group */

    fprintf(outf,"\nTransformed Ability distribution for New Group\n"
                 "Based on %s Method\n"
                 "Results stored in file named: %s\n\n",
                 (!strcmp(ST,"MS") ? "Mean/Sigma" :
                 (!strcmp(ST,"MM") ? "Mean/Mean" : 
                 (!strcmp(ST,"HA") ? "Haebara" : "Stocking-Lord"))),  
           DistNewFileST);

    inf = fopen(DistNewFileST,"r");
    fgets(line1,998,inf);                                 /* skip the first line */
    while(fgets(line1,998,inf) != NULL)
      fprintf(outf,"%s",strcat(line1,fgets(line2,998,inf))); 
    fclose(inf);
    fprintf(outf,"------------------------------------------------\n");
    
  }

  /* deallocate memory */

  StComItemDeAlloc(ComItems, &StInfo);
  StItemDeAlloc(NewItems, "new", &StInfo);
  StItemDeAlloc(OldItems, "old", &StInfo);
  StContDeAlloc(&StInfo);

  return;
}

/***********************************************************************************/

