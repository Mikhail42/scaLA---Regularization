package LAfunction.Vector

  object Error{
    /** norm(a.-b)/norm(a)
    */
    def error(a: Array[Double], b: Array[Double]): Double = {
      val c = BasicOperation.minus(a, b)
      Norm.norm1(c)/Norm.norm1(a)
    }
  }