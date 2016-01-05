package LAfunction.Vector

  
  object Inverse{
    /** rearranges a vector xs in accordance with the permutation vector q 
     *  @param x = (1,2)
     *  @param q = (1,0)
     *  @return (2,1)
     *  */
    def swapInverse(xs: Array[Double], q:Array[Int]): Unit = {
      for (i <- 0 until xs.length) {
        val ind = Search.indexOf(q,i)
        if (i!=ind) {
          val x = xs(i); xs(i)=xs(ind); xs(ind)=x;
          val y = q(i); q(i)=q(ind); q(ind)=y;
        }
      }        
    }
  }