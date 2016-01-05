package LAfunction.Vector
object BasicOperation{
    
    /** a(:).+b(:) */
    def plus(a: Array[Double], b: Array[Double]): Array[Double] = {
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)+b(i); c
    }
    
    /** a(:).-b(:) */
    def minus(a: Array[Double], b: Array[Double]): Array[Double] = {
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)-b(i); c
    }
    
    
    /** f*a(:) */
    def multi(a: Array[Double], f: Double):  Array[Double] = {
      import concurrent.Future
      import concurrent.ExecutionContext.Implicits.global
      val c = new Array[Double](a.length); for (i <- 0 until a.length) c(i)=a(i)*f;  c
    }
    
    /** a(i1:i2-1).*b(i1:i2-1) = sum(ai*bi : i in [i1,i2)) */
    def scalar(a: Array[Double], b: Array[Double], i1: Int, i2: Int): Double = {
      var sum = 0.0; for (i <- i1 until i2) sum+=a(i)*b(i); sum
    }
    
    /** a(i1:i2-1).*B(i1:i2-1)(j) */
    def scalar(a: Array[Double], B: Array[Array[Double]], j: Int, i1: Int, i2: Int): Double = {
      var sum = 0.0; for (i <- i1 until i2) sum+=a(i)*B(i)(j); sum
    }
    
    /** a(:).*a(:) */
    def product(a: Array[Double]): Double = {
      var mult = 1.0; for (x <- a) mult*=x; mult
    }
    
    /** 
     *  inner product
     *  multi(a, b,0,a.length) 
     **/ 
    def scalar(a: Array[Double], b: Array[Double]): Double = scalar(a, b,0,a.length)
    
    /** multi(a,B,j,0,a.length) */
    def scalar(a: Array[Double], B: Array[Array[Double]], j: Int): Double = scalar(a,B,j,0,a.length)
    
    def scalar(A: Array[Array[Double]], k: Int, B: Array[Array[Double]], l: Int): Double = {
      var sum = 0.0
      for (i <- 0 until A.length) sum += A(i)(k)*B(i)(l)
      sum
    }
    
    /**
     * outer product
     *  a*b^T 
     **/
    def multi(a: Array[Double], b: Array[Double]): Array[Array[Double]] = {
      val res = Array.ofDim[Double](a.length,b.length);
      for (i <- 0 until a.length) for (j <- 0 until b.length) res(i)(j) = a(i)*b(j)
      res
    }
  }  