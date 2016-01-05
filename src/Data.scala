
object Data {
  
  import BasicData._
  import Alpha._
  import Beta._
  import SLAE._
  import LAfunction.MatrixFunction.SomeOperation._
  import LAfunction.MatrixFunction.Multiply._
  import LAfunction.MatrixFunction.Info._
  import LAfunction.Vector.Search._
  import LAfunction.Number.DoubleFunction._
  import LAfunction._
  import math._
  
  /**
   * using constants
   */
  object BasicData{
    val in = new java.util.Scanner(System.in)
    println("Введите a, b, N")
    val a = in.nextDouble() //0.0
    val b = in.nextDouble()
    val N = in.nextInt
    val h = (b-a)/N.toDouble
    val size = N+1
    
    val name = ("Решение интегрального уравнения "
              +"\u222B"+"sec(ts)"+"\u03D5"+"(s)ds=4t, s:="
              +a+":"+h+":"+b)

    val xs = new Array[Double](size)
    xs(0) = a; for (i <- 1 until size) xs(i) = xs(i-1) + h
  
    val ts = new Array[Double](N+1)
    ts(0) = a; for (i <- 1 until size) ts(i) = ts(i-1) + h
    
    val K = Array.ofDim[Double](size, size)
    for (i <- 0 until size; j <- 0 to i)         K(i)(j) = sec(xs(i)*ts(j))
    for (i <- 0 until size; j <- i+1 until size) K(i)(j) = K(j)(i)
    
    def f(t: Double) = four(t)
  }
  
  /**
   * using BasicData
   */
  object As{
    
    lazy val as = {
      val res = new Array[Double](size)
      
      val oneDivH = 1.0/h
      val hDiv3   = h/3.0
      // Simpson's rule
      if (isDiv2(N)) {
        res(0) =  hDiv3; res(N) = hDiv3
        val fourHDiv3 = four(hDiv3)
        for (i <- 1 until N by 2) res(i) = fourHDiv3
        val twoHDiv3 = two(hDiv3)
        for (i <- 2 until N by 2) res(i) = twoHDiv3
      }
      // Trapezoidal rule
      else {
        val hDiv2 = h*0.5; res(0) = hDiv2; res(N) = hDiv2
        for (i <- 1 until N) res(i) = h
      }
      
      res
    }
  }
  
  /**
   * using constants
   */
  object Alpha{
    
    val nSqrAlphas = 20
    lazy val sqrAlphas = { 
      val res = new Array[Double](nSqrAlphas)
      for (i <- 0 until nSqrAlphas) res(i) = math.pow(10.0,-i)
      res
    }
    
    val nBestSqrAlphas = 4
    lazy val bestSqrAlphas = {
      val res = new Array[Double](nBestSqrAlphas)
      for (i <- 0 until nBestSqrAlphas) 
        res(i) = sqrAlphas((i<<1)+1)
      res
    }
    
  }
 
  /**
   * using Alpha, Eignvalues (Alpha)
   */
  object Beta{
    
    val namesBeta = {
      val res = new Array[String](3)
      res(0) = "1"
      res(1) = "\u03B1"
      res(2) = "sqrt(\u03B1^2+0.5*\u03C3_m^2(A))"
      res
    }
    
    
    val nBetas = 3     // don't edit!!! see down
    lazy val betas = {
      val res = new Array[Array[Double]](nBetas)
      
      res(0) = new Array[Double](1)
      res(0)(0) = 1.0
      
      res(1) = new Array[Double](nBestSqrAlphas)
      for (j <- 0 until nBestSqrAlphas) res(1)(j) = sqrt(bestSqrAlphas(j))
      
      res(2) = new Array[Double](nBestSqrAlphas)
      for (j <- 0 until nBestSqrAlphas) {
        //val value = getMaxAndMin(Eignvalues.matEigns(j))._2
        res(2)(j) = sqrt(sigma_m*0.5+Alpha.bestSqrAlphas(j))
      }
      
      res
    }
    
  }
  
  /**
   * using Alpha, Cond (Alpha)
   */
  object Eignvalues{
     
    /**
     * using Alpha
     */
    lazy val matEigns = {
      val res = new Array[Array[Double]](Alpha.nBestSqrAlphas)
      for (j <- 0 until nBestSqrAlphas) { 
        try {
          val eigns = eigenvaluesA_k(ATAsWithAlpha(j))
          res(j) = eigns._1
          val exception = eigns._2
        } catch {
          case e: Exception => {}
        }
      }
      res
    }
    
    /**
     * using Alpha, Eigenvalues (Alpha), Cond (Alpha), SLAE (BasicData)
     */
    lazy val superMatEigns = {      
      val res = new Array[Array[Array[Double]]](3)
      
      for (i <- 0 until 3) res(i) = new Array[Array[Double]](Alpha.nBestSqrAlphas)
      
      def getEigns(alpha: Double, beta: Double): Array[Double] = {
        val elem = -alpha/beta
        
        for (l <- 0 until size)        SLAE.resultMat(l)(l) += beta
        for (l <- size until twoSize)  SLAE.resultMat(l)(l) += elem
        
        val res = eigenvaluesA_k(SLAE.resultMat)._1
        
        for (l <- 0 until size)        SLAE.resultMat(l)(l) -= beta
        for (l <- size until twoSize)  SLAE.resultMat(l)(l) -= elem
        
        res
      }
      
      for (j <- 0 until  Alpha.nBestSqrAlphas){
        res(0)(j) = getEigns(Alpha.bestSqrAlphas(j), betas(0)(0))
        res(1)(j) = getEigns(Alpha.bestSqrAlphas(j), betas(1)(j))
        res(2)(j) = getEigns(Alpha.bestSqrAlphas(j), betas(2)(j))
      }
      
      res 
    }
  }
  
  /**
   * using Alpha, SLAE (BasicData)
   */
  object Cond{
    
    def getCond(eigns: Array[Double]): Double = {
      val maxMin = getMaxAndMin(eigns)
      maxMin._1/maxMin._2
    }
    
    lazy val condA = cond2(A)
    
    lazy val condATA = cond2(ATA)
    
    
    /**
     * using Alpha
     */
    lazy val vectCond = { 
      val res = new Array[Double](nBestSqrAlphas)
      for (i <- 0 until nBestSqrAlphas) res(i) = getCond(Eignvalues.matEigns(i))
      res
    }
    
    /**
     * using Alpha, Eignvalues (Alpha), Beta (constants), Cond (Alpha), SLAE (BasicData)
     */
    lazy val matCond = {
      val res = Array.ofDim[Double](Beta.nBetas, Alpha.nBestSqrAlphas)
      for (i <- 0 until Beta.nBetas; j <- 0 until  Alpha.nBestSqrAlphas)
        res(i)(j) = getCond(Eignvalues.superMatEigns(i)(j))
      res
    }
  
  }
  
  /**
   * using BasicData, As 
   */
  object SLAE{
    
    import As._
    
    /**
     * A*x = b
     */
    lazy val A = MatrixFunction.Multiply.scalarMulti(K, as) 
      
    private val n = A.length
    /*private lazy val x = {
      val res = new Array[Double](n)
      for (i <- 0 until n) res(i) = i+1
      res
    }*/ 
      // multi(A, x)
    lazy val b = {
      val res = new Array[Double](n)
      for (i <- 0 until n) res(i) = f(xs(i))
      res
    }
    
    /**
     * for 
     * 	(A^T A + (alpha^2)*E)x = A^T b
     */
    lazy val ATA = multiATA(SLAE.A)
    lazy val ATb = multiATb(SLAE.A, SLAE.b)
    lazy val sigma_m = getMaxAndMin(eigenvaluesA_k(SLAE.resultMat)._1)._2
    
    /**
     * for
     * 	(b   A     ) (r/b) = (b)
     * 	(AT  -a^2/b) (x  ) = (0)
     */
    val twoSize = (size<<1)
    lazy val resultMat = {
      val res = Array.ofDim[Double](twoSize, twoSize)
      for (k <- 0 until size; l <- size until twoSize) res(k)(l) = SLAE.A(k)(l-size)
      for (k <- size until twoSize; l <- 0 until size) res(k)(l) = SLAE.A(l)(k-size)
      res
    }
    lazy val ATAsWithAlpha = {
      val res = new Array[Array[Array[Double]]](Alpha.nBestSqrAlphas)
      for (i <- 0 until Alpha.nBestSqrAlphas){
        res(i) = Array.ofDim[Double](n,n)
        LAfunction.MatrixFunction.Copy.copy(ATA, res(i))
        plusDiag(res(i), bestSqrAlphas(i))
      }
      res
    }
    lazy val rightPart = {
      val res = new Array[Double](twoSize)
      for (i <- 0 until size) res(i) = SLAE.b(i)
      res
    }
    
    lazy val ATAsWithAlphaAndBeta = {
      
      val res = new Array[Array[Array[Double]]](Beta.nBetas)
      for (i <- 0 until Beta.nBetas){
        res(i) = Array.ofDim[Double](twoSize, twoSize)
        LAfunction.MatrixFunction.Copy.copy(resultMat, res(i))
        val alpha = Alpha.bestSqrAlphas(Solution.betstIndAlpha)
        val betaElem = Beta.betas(i)(if (i==0) 0 else Solution.betstIndAlpha)
        val eleme = -alpha/betaElem
        for (l <- 0 until twoSize)       res(i)(l)(l) += betaElem
        for (l <- twoSize until twoSize) res(i)(l)(l) += eleme
      }
      res
    }
  }
  
  object Solution {
    
    import LAfunction.MatrixFunction.Solution._
    
    val betstIndAlpha = 3 // 0 <= betstIndAlpha < nBestAlphas
    
    object Jacobi{
      lazy val solveATAsWithAlpha = {
        val res = Array.ofDim[Double](Alpha.nBestSqrAlphas, ATA.length)
        for (i <- 0 until Alpha.nBestSqrAlphas) {
          val method = new LAfunction.JacobiAndSeidel(ATAsWithAlpha(i), ATb)
          res(i) = method.xsSeidel 
        }
        res
      }
      
      lazy val solveATAsWithAlphaAndBeta = {
        val res = Array.ofDim[Double](Beta.nBetas, xs.length)
        for (i <- 0 until Beta.nBetas){
          val method = new LAfunction.JacobiAndSeidel(ATAsWithAlphaAndBeta(i), rightPart)
          Array.copy(method.xsSeidel, xs.length-1, res(i), 0, xs.length)
        }
        res
      }
    }
    
    object LU {
      
      lazy val solveATAsWithAlpha = {
        val res = Array.ofDim[Double](Alpha.nBestSqrAlphas, ATA.length)
        for (i <- 0 until Alpha.nBestSqrAlphas) 
          res(i) = LAfunction.MatrixFunction.Solution.solveLU(ATAsWithAlpha(i), ATb)
        res
      }
    
      lazy val solveATAsWithAlphaAndBeta = {
        val res = Array.ofDim[Double](Beta.nBetas, xs.length)
        for (i <- 0 until Beta.nBetas){
          val solve = solveLU(ATAsWithAlphaAndBeta(i), rightPart)
          Array.copy(solve, xs.length-1, res(i), 0, xs.length)
        }
        res
      }
      
    }
    
    object LUP {
      
      lazy val solveATAsWithAlpha = {
        val res = Array.ofDim[Double](Alpha.nBestSqrAlphas, ATA.length)
        for (i <- 0 until Alpha.nBestSqrAlphas) 
          res(i) = LAfunction.MatrixFunction.Solution.solveLUP(ATAsWithAlpha(i), ATb)
        res
      }
    
      lazy val solveATAsWithAlphaAndBeta = {
        val res = Array.ofDim[Double](Beta.nBetas, xs.length)
        for (i <- 0 until Beta.nBetas){
          val solve = solveLUP(ATAsWithAlphaAndBeta(i), rightPart)
          Array.copy(solve, xs.length-1, res(i), 0, xs.length)
        }
        res
      }
      
    }
    
    object QR {
    
      lazy val solveATAsWithAlpha = {
        val res = Array.ofDim[Double](Alpha.nBestSqrAlphas, ATA.length)
        for (i <- 0 until Alpha.nBestSqrAlphas) 
          res(i) = solveQR(ATAsWithAlpha(i), ATb)
        res
      }
      
      lazy val solveATAsWithAlphaAndBeta = {
        import LAfunction.MatrixFunction.Solution._
        val res = Array.ofDim[Double](Beta.nBetas, xs.length)
        for (i <- 0 until Beta.nBetas){
          val solve = solveQR(ATAsWithAlphaAndBeta(i), rightPart)
          Array.copy(solve, xs.length-1, res(i), 0, xs.length)
        }
        res
      }
    }
  }
}