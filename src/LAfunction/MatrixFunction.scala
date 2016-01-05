package LAfunction

import math._
import LAfunction.Vector.BasicOperation._
import LAfunction.Number.DoubleFunction

/**
 * @author Ionkin Mikhail

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
object MatrixFunction {

  object Multiply {
    
    /**
     * A^T A
     */
    def multiATA(A: Array[Array[Double]]): Array[Array[Double]] = {
      val n = A.length
      val res = Array.ofDim[Double](n,n)
      for (i <- 0 until n; j <- 0 until n)
        res(i)(j) = LAfunction.Vector.BasicOperation.scalar(A,i, A,j) 
      res
    }
    
    def multiATb(A: Array[Array[Double]], b: Array[Double]): Array[Double] = {
      val n = A.length
      val res = new Array[Double](n)
      for (j <- 0 until n) res(j) = LAfunction.Vector.BasicOperation.scalar(b, A, j)
      res
    }
    
    def multi(A: Array[Array[Double]], f: Double, size: Int): Array[Array[Double]] = {
      val n = size
      val res = Array.ofDim[Double](n,n)
      for (i <- 0 until n; j <- 0 until n) res(i)(j) = A(i)(j)*f
      res
    }
    
    def multi(A: Array[Array[Double]], f: Double): Array[Array[Double]] = multi(A,f,A.length)
    
    def multi(A: Array[Array[Double]], b: Array[Double], size: Int) : Array[Double] = {
      val c = new Array[Double](size)
      for(i <- 0 until size) c(i) = LAfunction.Vector.BasicOperation.scalar(A(i),b,0,size)
      c
    }
    
    def scalarMulti(A: Array[Array[Double]], b: Array[Double]): Array[Array[Double]] = {
      val res = Array.ofDim[Double](A.length,A.length)
      for (i<- 0 until A.length; j <-0 until A.length) res(i)(j) = A(i)(j)*b(j)
      res
    }
    
    def multi(A: Array[Array[Double]], b: Array[Double]) : Array[Double] = multi(A,b,A.length)
     
    /**
     * size(A) = (n,n) = size(B), size in \doubleZ
     */
    def multi(A: Array[Array[Double]], B: Array[Array[Double]], size: Int): Array[Array[Double]] = {
      val n = size
      val C =  Array.ofDim[Double](n,n)
      for (i <- 0 until n) for (j <- 0 until n) C(i)(j) = LAfunction.Vector.BasicOperation.scalar(A(i), B, j, 0, n)
      C
    }   
    
    /**
     * size(A) = (n,n) = size(B)
     * multi(A,B,A.length)
     */
    def multi(A: Array[Array[Double]], B: Array[Array[Double]]): Array[Array[Double]] =
        multi(A,B,A.length);
    
  }
  
  object EditMatrix{
    
    def swapColumn(mat: Array[Array[Double]], i1: Int, i2: Int): Unit = {
      if (i1 != i2) {
        for (i <- 0 until mat.length){ 
          mat(i)(i1) += mat(i)(i2) 
          mat(i)(i2) = mat(i)(i1)-mat(i)(i2)
          mat(i)(i1) -= mat(i)(i2)
        }
      }
    }
    
    def swapString(mat: Array[Array[Double]], i1: Int, i2: Int): Unit = {
      if (i1 != i2) {
        for (j <- 0 until mat.length){ 
          val x = mat(i1)(j)
          mat(i1)(j) = mat(i2)(j)
          mat(i2)(j) = x
        }
      }
    }
    
    def transp(mat: Array[Array[Double]]): Unit = {
      val n = mat.length
      for (i <- 0 until n)
            for (j<- 0 until i) {
              val x = mat(i)(j);
              mat(i)(j) = mat(j)(i)
              mat(j)(i) = x
            }
    }
    
  }
 
  object Info{
    
    def detLUviaU(U: Array[Array[Double]]) = SomeOperation.detDiagMatr(U)
    
    def det(A:Array[Array[Double]]) = {
      val LU = Decomposition.getLU(A, false, null)
      detLUviaU(LU._2)
    }
    
    def cond(A: Array[Array[Double]], L: Array[Array[Double]], U: Array[Array[Double]]) = 
        normMS(A)*normMS(Inverse.inverseMatrix(L, U))
    
    def cond(A: Array[Array[Double]]) = normMS(A)*normMS(Inverse.inverseMatrix(A)) 
    
    /**
     * p-norm, where p = infinity
     * the maximum absolute row sum of the matrix
     */
    def normMS(A: Array[Array[Double]], size: Int): Double = {
      import LAfunction.Vector.Norm.norm1
      var maxA = norm1(A(0),size); 
      for (i <- 1 until size) maxA = max(maxA,norm1(A(i),size)); maxA
    }
    
    def normMS(a: Array[Array[Double]]): Double = normMS(a,a.length);
    
    /**
     * p-norm, where p = 1
     * the maximum absolute column sum of the matrix
     */
    def normMC(A: Array[Array[Double]], size: Int): Double = {
      import LAfunction.Vector.Norm.norm1
      var maxA = norm1(A,0,size); for (i <- 1 until size) maxA = max(maxA,norm1(A,i,size)); maxA
    }
    
    def normMC(A: Array[Array[Double]]): Double = normMC(A,A.length);
    
    /**
     * (p-norm, where p = 2)
     */
    def normSpectr(eigenvalues: Array[Double]): Double = {
      val maxmin = LAfunction.Vector.Search.getMaxAndMin(eigenvalues)
      sqrt(maxmin._1/maxmin._2)
    }
    
    
    var QR_ERROR = 1e-4
    var QR_ITER  = 1000
    
    private def getA_k(A: Array[Array[Double]]): (Array[Array[Double]], java.lang.IllegalArgumentException)  = {  
      val n = A.length      
      var QR = Decomposition.QR.getQR_HR(A)
      var curA = Multiply.multi(QR._2, QR._1); // for (i <- 0 until n) curA(i)(i)+=Ann
      var norm2 = normMS(curA)
      var norm1 = norm2+QR_ERROR*10+1.0
      var iter = QR_ITER
      while (abs(norm2-norm1)>QR_ERROR && iter>0){
        // var Ann = updA(n-1)(n-1); for (i <- 0 until n) updA(i)(i)-=Ann
        QR = Decomposition.QR.getQR_GR(curA)
        curA = Multiply.multi(QR._2, QR._1); // for (i <- 0 until n) curA(i)(i)+=Ann
        norm1 = norm2
        norm2 = normMS(curA)
        iter-=1
      }
      if (iter == 0) 
        return(curA, new java.lang.IllegalArgumentException(
          "Последовательность {Ak}, k=0:infinity не сходится по норме на "+QR_ITER+" итерациях (с точностью до QR_ERROR = "+QR_ERROR+").\n"+
          "Предполагаемые собственные числа: "+
              LAfunction.Vector.Print.toString(
                {
                  val diag = new Array[Double](n)
                  for (i <- 0 until n) diag(i) = curA(i)(i)
                  diag
                })))
      else 
        (curA, null)
    }
    
    def eigenvaluesA_k(A: Array[Array[Double]]): (Array[Double], java.lang.IllegalArgumentException) = {  
      val n = A.length
      val res = getA_k(A: Array[Array[Double]]) 
      val Ak = res._1
      val exception = res._2
      
      val diag = new Array[Double](n)
      for (i <- 0 until n) diag(i) = Ak(i)(i)
      (diag, exception)
    }
          
    def cond2(A: Array[Array[Double]]): (Double, java.lang.IllegalArgumentException) = {
      val res = eigenvaluesA_k(A)
      val resEign = res._1
      val exception = res._2
        
      val maxAndMin = LAfunction.Vector.Search.getMaxAndMin(resEign)
      (maxAndMin._1/maxAndMin._2, exception)
    } 

  }
  
  object Inverse {
      
    def inverseTriangleUp(U: Array[Array[Double]]): Array[Array[Double]] = {
        val n = U.length
        val res = Array.ofDim[Double](n,n)
        //[http://alexandr4784.narod.ru/kaplan5/kaplan_5_06.pdf p. 146] 
        for (i <- 0 until n) res(i)(i) = 1.0/U(i)(i);
        for (i <- 0 until n)  for (j <- i+1 until n) 
          res(i)(j) = -res(j)(j)*LAfunction.Vector.BasicOperation.scalar(res(i), U, j, 0, j)
        res
    }
    
    def inverseTriangleDown(L: Array[Array[Double]]): Array[Array[Double]] = {
        val n = L.length
        val res = Array.ofDim[Double](n,n)
        //[http://alexandr4784.narod.ru/kaplan5/kaplan_5_06.pdf p. 149] 
        for (i <- 0 until n) res(i)(i) = 1.0/L(i)(i);
        for (i <- 0 until n) for (j <- 0 until i) 
          res(i)(j) = -res(i)(i)*LAfunction.Vector.BasicOperation.scalar(L(i),res, j, 0, i)
        res
    }
    
    /**
    * A = LU; 
    * [1, eqs. 1.8-1.10]
    * @return A^(-1)
    */
    def inverseMatrix(L: Array[Array[Double]], U: Array[Array[Double]]): Array[Array[Double]] = {
      val n = L.length
      import LAfunction.Vector.BasicOperation._
      val mat = Array.ofDim[Double](n,n);
      /** [1, eqs. 1.8-1.10] */
      def setMat(i: Int, j: Int) : Unit = {
        mat(i)(j) = 
          if (i == j)    (1.0 - scalar(U(i), mat,j,j+1,n))/U(j)(j)
          else if (i < j)    (- scalar(U(i), mat,j,i+1,n) /U(i)(i))
          else               (- scalar(U(i), L,  j,j+1,n))
      }
      
      var i = n-1
      while (i>=0){
        var j = n-1
        while (j>=0) {
          setMat(i,j); j-=1
          }
        i-=1
      }
      
      mat
    }
    
    def inverseMatrix(A: Array[Array[Double]]):  Array[Array[Double]] = {
      val LU = Decomposition.getLU(A, false, null)
      inverseMatrix(LU._1, LU._2)
    }
  }
  
  object Solution{
    
    /** 
     *  solution of SLAE Ax = b 
     *  не тестировалась!!!
     *  */
    def solveLUP(A: Array[Array[Double]], b: Array[Double]): Array[Double] =  { 
      val q = new Array[Int](A.length)
      for (i <- 0 until A.length) q(i) = i
      val LU = Decomposition.getLU(A, true, q)
      val x = Solution.solveLU(LU._1, LU._2, b);
      LAfunction.Vector.Inverse.swapInverse(x,q); 
      x
    }
    
    /**
    * [1], p. 7, eq. (1.5)
    */
    def solveLU(L: Array[Array[Double]], U: Array[Array[Double]], b: Array[Double]): Array[Double] = {
      val n = L.length
      
      val ys = new Array[Double](n)
      for (i <- 0 until n) ys(i) = b(i) - scalar(L(i),ys,0,i)
      
      val xs = new Array[Double](n)
      var i = n-1
      while (i>=0) {
        xs(i) = (1.0/U(i)(i))*(ys(i) - scalar(U(i),xs,i+1,n))
        i-=1
      }
      
      xs
    }
    
    /**
    * [1], p. 7, eq. (1.5)
    */
    def solveLU(A: Array[Array[Double]], b: Array[Double]): Array[Double] = {
      val LU = Decomposition.getLU(A, false, null)
      solveLU(LU._1, LU._2, b)
    }
    
    def solveQR(QR: (Array[Array[Double]], Array[Array[Double]]), b: Array[Double]): Array[Double] = {
      val n = QR._1.length
      val Q = QR._1
      val R = QR._2
     
      val newB = new Array[Double](n)
      for (k <- 0 until n) newB(k) = scalar(b,Q,k)
      
      val xs = new Array[Double](n)
      for (i <- n-1 to 0 by -1)
        xs(i) = (1.0/R(i)(i))*(newB(i) - scalar(R(i),xs,i+1,n))      
      xs
    }
    
    def solveQR(A: Array[Array[Double]], b: Array[Double]): Array[Double] = 
      solveQR(Decomposition.QR.getQR_HR(A), b)
     
  }
   
  object Decomposition{
  
    object QR{
      
      /** 
      *  Householder reflection 
      *  [1: 1.3--Метод отражений]
      *  */
      def getQR_HR(A: Array[Array[Double]]): (Array[Array[Double]], Array[Array[Double]]) = {
        val n = A.length
           
        var curH1: Array[Array[Double]] = null; 
        val curH2 = Array.ofDim[Double](n,n)
        var curR = Array.ofDim[Double](n,n); 
        Copy.copy(A, curR)
          
        solve()
           
         def solve(){
          for (k <- 0 until n-1){              
            val v = new Array[Double](n-k); 
            for (i <- k until n) v(i-k) = curR(i)(k); 
            // p = v; [1, 1.3 -- eq 1.14]
            v(0) += DoubleFunction.specSign(v(0))*LAfunction.Vector.Norm.norm1(v)   
                
            curH1 = getHouseholdersMatrix(v)
            setH2()
    
            curR = Multiply.multi(curH2,curR)            
          }
          
          /**
          * Householders Matrix
          */
          def getHouseholdersMatrix(p: Array[Double]): Array[Array[Double]] = {
            val matrix = Multiply.multi(multi(p,p), -2.0/scalar(p,p))
            for (i <- 0 until matrix.length) matrix(i)(i) += 1
            matrix
          }
          
          def setH2(){
            val i2 = n - curH1.length
            for (i <- 0 until i2; j <- 0 until i2) 
              curH2(i)(j) = if (i == j) 1.0 else 0.0 
            for (i <- 0 until i2; j <- i2 until n) curH2(i)(j) = 0.0
            for (j <- 0 until i2; i <- i2 until n) curH2(i)(j) = 0.0
            for (i <- i2 until n; j <- i2 until n) curH2(i)(j) = curH1(i-i2)(j-i2)
          }
          
        }
        val Q = Multiply.multi(A, Inverse.inverseTriangleUp(curR));
        (Q, curR)
      }
      
      // Givens rotations
      def getQR_GR(A: Array[Array[Double]] ): (Array[Array[Double]],Array[Array[Double]]) = { 
        val n = A.length
        var curR = Array.ofDim[Double](n,n)
        Copy.copy(A, curR)        
        
        for (i <- 0 until n){
          for (j <- 0 until i) 
            if (!DoubleFunction.isZero(curR(i)(j))) {
              val colm = new Array[Double](n); for (k <- 0 until n) colm(k)=curR(k)(j)
              val (c,s) = getGivensMatrixCS(colm,i,j)
              curR = multiGR(curR,c,s,i,j)
            }
        }
        
        /** 
        *  (c,s)
        **/
        def getGivensMatrixCS(v: Array[Double], i: Int, j: Int): (Double, Double) = {
          val n = v.length
          val m1tay = 1.0/sqrt(v(i)*v(i)+v(j)*v(j))
          (v(j)*m1tay, v(i)*m1tay)
        }
        
        /**
        * 1 0  0 0    11 12 13 14   11      12      13      14 
        * 0 c  s 0  * 21 22 23 24 = c21+s31 c22+s32 c23+s33 c24+s34 
        * 0 -s c 0    31 32 33 34  -s21+c31-s22+c32
        * 0 0  0 1    41 42 43 44   41      42      43      44
        **/
        def multiGR(B: Array[Array[Double]], c: Double, s: Double,  i00: Int, j00: Int): Array[Array[Double]] = {
          val n = A.length
          val res = Array.ofDim[Double](n,n)
          for (i <- 0 until n) if (i!=i00 && i!=j00) for (j <- 0 until n) res(i)(j) = B(i)(j)
          var i0 = min(i00,j00); var j0 = max(i00,j00); 
          for (j <- 0 until n) res(i0)(j) =  c*B(i0)(j)+s*B(j0)(j)
          for (j <- 0 until n) res(j0)(j) = -s*B(i0)(j)+c*B(j0)(j)
          res
        }
        
        val Q = Multiply.multi(A, Inverse.inverseTriangleUp(curR))
        (Q, curR)
      }
  
    }
    
    def getLLT(matr: Array[Array[Double]]): Array[Array[Double]]= {
      val n = matr.length
      import math._
      val L = Array.ofDim[Double](n,n)
      for (i <- 0 until n)
        for (j <- 0 to i)
          if (i==j) L(i)(i) = sqrt(matr(i)(i)- scalar(L(i),L(i),0,i))
          else L(i)(j) = (1.0/L(j)(j))*(matr(i)(j) - scalar(L(i),L(j),0,j))
      L
    }
    
    /**
     * LU-decomposition. 
     * [1 : 1.1,1.2], [2: 1.2, 4.4] 
     * @param dwpFlag -- (true => decomposition with pivoting, false => simple decomposition)
     * @q -- the resulting vector permutations (if dwpFlag == true). 
     * @return (L,U). q is UPDATE
     */
    @throws(classOf[java.lang.IllegalArgumentException])
    def getLU(matr: Array[Array[Double]], dwpFlag: Boolean, q: Array[Int]):  (Array[Array[Double]], Array[Array[Double]]) = {       
      val n = matr.length
      val U = Array.ofDim[Double](n,n)
      val L = Array.ofDim[Double](n,n)
      /**
       *  [2]: 
       *   стр. 11: Можно чередовать вычисление строк U(k)(:) и столбцов L(:)(k)
       *   стр. 42: Устойчивость не гарантирована 
       */
    /*  for (k <- 0 until n) {
        for (j <- k until n) U(k)(j) = matr(k)(j) - scalar(L(k),U,j,0,k);
        if (dwpFlag) {
          val ind = LAfunction.Vector.Search.indexOfMaxAbs(U,k,k) 
          EditMatrix.swapString(matr, k, ind)
          EditMatrix.swapString(U, k, ind)
          if (k!=ind) { 
            val g = q(k)
            q(k) = q(ind)
            q(ind) = g
          }
        }
      */
      for (k <- 0 until n) {
        for (j <- k until n) U(k)(j) = matr(k)(j) - scalar(L(k),U,j,0,k);
        if (dwpFlag) {
          // индекс элемента с наибольшим модулем
          // !!!!!!!!!!!!!!!!!!!!!!!!!!!
          // val ind = LAfunction.Vector.indexOfMaxAbs(matr(k),k) // val ind = LAfunction.Vector.indexOfMaxAbs(U(k),k)
          val ind = LAfunction.Vector.Search.indexOfMaxAbs(matr(k),k) 
          EditMatrix.swapColumn(matr, k, ind)
          EditMatrix.swapColumn(U, k, ind)
          if (k!=ind) { 
            val g = q(k)
            q(k) = q(ind)
            q(ind) = g
          }
        }
        if (DoubleFunction.isZero(U(k)(k))) throw new java.lang.IllegalArgumentException(
              "U("+k+")("+k+")="+U(k)(k)
              +". \n\t\t\t LU-разложение невозможно ввиду необходимости деления на сверхмалое число. \n "
              +"Вектор перестановок: " + LAfunction.Vector.Print.toString(q)+'\n'
              +k+"-я строка: " + LAfunction.Vector.Print.toStringWithE(U(k))
            )
        val m1Ukk = 1.0/U(k)(k)
        L(k)(k) = 1.0
        for (i <- k+1 until n) L(i)(k) = (matr(i)(k) - scalar(L(i),U,k,0,k))*m1Ukk
      }      
      (L,U)
    }
  }
  
  object SomeOperation{
    
    def plusDiag(A: Array[Array[Double]], x: Double) = 
      for (i <- 0 until A.length) A(i)(i) += x
      
    def detDiagMatr(A: Array[Array[Double]]): Double = {
      var res = 1.0; for (i <- 0 until A.length) res*=A(i)(i); res 
    }
    
    def countZeros(A: Array[Array[Double]]): Int = {
      var count = 0; for (vect <- A; elem <- vect) 
        if (DoubleFunction.isZero(elem)) count+=1; count;
    } 
  }
  
  object IO{
        
    def toString(A: Array[Array[Double]]): String = {
      var maxLenght = 0
      val format = "%."+2+'f'
      for (str <- A; elem <- str) maxLenght = max(maxLenght, format.format(elem).length)
      toString(A, maxLenght, 2)
    }
    def toString(A: Array[Array[Double]], count1: Int, count2: Int): String = {
      val out = new java.lang.StringBuilder()
      for (i <- 0 until A.length) 
        out.append(LAfunction.Vector.Print.toString(A(i), count1, count2)).append('\n')
      out.toString()
    }
    
    def input(strs: Array[String]): Array[Array[Double]] = {
      val in1 = strs(0).split(' ').map{ x => x.toDouble }
      val n = in1.length
      val res =  Array.ofDim[Double](n,n);  
      Array.copy(in1,0,res(0),0,n);
      for (i <- 1 until n) res(i) = strs(i).split(' ').map{ x => x.toDouble }
      res
    }
  }
  
  object Copy {
    
    def copy(inp: Array[Array[Double]], out: Array[Array[Double]], size: Int): Unit = {
      for (i <- 0 until size) Array.copy(inp(i), 0, out(i), 0, size)
    }
    
    def copy(inp: Array[Array[Double]], out: Array[Array[Double]]): Unit = copy(inp,out,inp.length)
  }
}