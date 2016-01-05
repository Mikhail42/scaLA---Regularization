package LAfunction.Vector

 object Print { 
    
    /** 
     *  @count the number of digits after the decimal point. 
     *  format = '%'+@count.toString+'f'.  
     *  @return sum(format.format(a(:))+' ')
     *  */
    def toString(a: Array[Double], count: Int): String = {
      val format = "%."+count+'f'
      val out = new java.lang.StringBuilder()
      for (i <- 0 until a.length) out.append(format.format(a(i))).append(' ')
      out.toString()
    }
    
    /** 
     *  @count the number of digits after the decimal point. 
     *  format = '%'+@count.toString+'f'.  
     *  @return sum(format.format(a(:))+' ')
     *  */
    def toString(a: Array[Double], count1: Int, count2: Int): String = {
      val out = new java.lang.StringBuilder()
      val format = "%"+count1+"."+count2+"f" 
      for (i <- 0 until a.length-1) out.append(format.format(a(i))).append(' ')
      out.append(format.format(a(a.length-1)))
      out.toString()
    }
    
    def toString(a: Array[Int]): String = {
      val out = new java.lang.StringBuilder()
      for (i <- 0 until a.length) out.append(a(i)).append(' ')
      out.toString()
    }
    
    def toString(a: Array[Double]): String = toString(a, 2)
    
    def toStringWithE(a: Array[Double]): String = {
      val out = new java.lang.StringBuilder()
      for (i <- 0 until a.length-1) out.append(a(i)).append(' ')
      out.append(a(a.length-1))
      out.toString()
    }
    
  }