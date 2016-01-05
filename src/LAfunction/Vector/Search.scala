package LAfunction.Vector

import math._

  object Search{
    
    /** 
     *  Search index maximal element in the abs(a(:)). Search begins with the beginIndex. For example: 
     *  @beginIndex = 1
     *  @param a (-11, -10, 9, 2)
     *  @return 1
     *  */ 
    def indexOfMaxAbs(a: Array[Double], beginIndex: Int): Int = {
      var curInd = beginIndex;
      for (i <- curInd+1 until a.length) if (abs(a(curInd)) < abs(a(i))) curInd = i; curInd    
    }  
    
    def indexOfMaxAbs(A: Array[Array[Double]], j: Int, beginIndex: Int): Int = {
      var curInd = beginIndex;
      for (i <- curInd+1 until A.length) 
        if (abs(A(i)(j)) < abs(A(curInd)(j))) curInd = i
      curInd    
    }  
    
    /** 
     *  @return (max(abs(a)),min(abs(a)))
     *  */
    def getMaxAndMin(a: Array[Double]): (Double,Double) = {
      var maxEV = abs(a(0)); var minEV = abs(a(0))
      for (i <- 1 until a.length){
        maxEV = max(maxEV,abs(a(i)))
        minEV = min(minEV,abs(a(i)))
      }
      (maxEV,minEV)
    }
    /**
     * It returns the index of the first matched element
     * @param a = (2, 3, 3, 4)
     * @param elem = 3
     * @return 1 
     */
    def indexOf[T: scala.reflect.ClassTag](a: Array[T], elem: T): Int = {
      for (i <- 0 until a.length) if (a(i).equals(elem)) return i
      -1
    } 
  }
  