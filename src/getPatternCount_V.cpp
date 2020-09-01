/*This function was fixed by Vlad Yatskou
30.06.2019
*/
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP getPatternCount_V (SEXP string, SEXP pattern)
{
  Rcpp::CharacterVector strVect(string);
  Rcpp::CharacterVector patVect(pattern);
  int rowsNum = strVect.size();
  int patternSize = strlen(patVect[0]);
  Rcpp:NumericVector result(rowsNum);
  
  for (int i = 0; i < rowsNum; i ++)
  {
    int stringSize = strlen(strVect[i]);
    result[i] = 0;
    
    bool matched = false;
    
    for (int stringI = 0; stringI + patternSize <= stringSize; stringI ++)
    {
      for (int patternJ = 0; patternJ < patternSize; patternJ ++)
      {
        matched = true;
        if (strVect[i][stringI + patternJ] != patVect[0][patternJ])
        {
          matched = false;
          break;
        }
      }
      if (matched)
      {
        result[i] ++;
		stringI += patternSize - 1;
      }
    }
  }

  return result;
}