options(scipen = 20)

library(jsonlite)
library(psych)
library(stringr)

originalExamplesPath = "~/Developer/EquatingRecipes/EqRCpp/tests/examples/resources/json/original_results/"
cppExamplesPath = "~/Developer/EquatingRecipes/EqRCpp/build/tests/examples/"

importCh2 = function(filename) {}

importCh3 = function(filename) {}

compareCh4 = function() {
  originalResults = fromJSON(paste0(originalExamplesPath, "chapter4.json"))
  
  cppResults = fromJSON(paste0(cppExamplesPath, "chapter4.json"))
  
  meanDiff = cppResults$EquatedRawScoreResults$equatedRawScoreMoments[, 1] -
    originalResults$equatedRawScoreMoments$Mean
}





mondatx = fromJSON("~/Developer/EquatingRecipes/EqRCpp/tests/examples/resources/json/datasets/mondatx.dat.json")
rawScores = mondatx[, (ncol(mondatx) - 1):ncol(mondatx)]



readUSTATS = function(linesRead) {
  startRow = 5
  # fprintf(fp,"\n\n%s\n\n",tt);
  # fprintf(fp,"Input filename:  %s\n\n",s->fname);
  # 
  # fprintf(fp,"Variable identifier: %c\n\n",s->id);
  # 
  # fprintf(fp,"Number of persons = %6d\n\n",s->n);
  # fprintf(fp,"Min score in data = %12.5f\n",s->mind);
  # fprintf(fp,"Max score in data = %12.5f\n",s->maxd);
  # fprintf(fp,"             Mean = %12.5f\n",s->mts[0]);
  # fprintf(fp,"             S.D. = %12.5f\n",s->mts[1]);
  # fprintf(fp,"             Skew = %12.5f\n",s->mts[2]);
  # fprintf(fp,"             Kurt = %12.5f\n\n",s->mts[3]);
  # 
  # fprintf(fp,"Min score for fd[] = %10.5f\n",s->min);
  # fprintf(fp,"Max score for fd[] = %10.5f\n",s->max); 
  # fprintf(fp,"Increment between scores = %10.5f\n",s->inc);
  # fprintf(fp,"Number of scores = %6d\n\n",s->ns); 
  # fprintf(fp,"       Score  Freq  CFreq    Rel Freq   Cum RFreq   Perc Rank\n\n");
  # for(i=0;i<=s->ns-1;i++) 
  #   fprintf(fp,"%12.5f%6d%7d%12.5f%12.5f%12.5f\n",score(i,s->min,s->inc),
  #           s->fd[i],s->cfd[i],s->rfd[i],s->crfd[i],s->prd[i]);
  # fprintf(fp,"\n");
  # for(i=1;i<=61;i++) fprintf(fp,"*");
}