package LAfunction

import MatrixFunction._
import MatrixFunction.Copy._
import MatrixFunction.Multiply._
import MatrixFunction.SomeOperation._
import LAfunction.Vector.BasicOperation._

/**
 * @author Ionkin Mikhail
 * This object implements the methods of Jacobi and Seidel
 * [1: 3.2,3.3]

 * [1] -- "С. Ю. Гоголева, Матрицы и вычисления, 2004
 * [2] -- "В. Б. Андреев, Численные методы"
 * [3] -- "J. E. Roman, The QR Algorithm. http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter3.pdf"
 * [4] -- "G. Fasshauer  http://www.math.iit.edu/~fass/477577_Chapter_11.pdf"
 */
class JacobiAndSeidel(A: Array[Array[Double]], b: Array[Double]) {
  
      private val n = A.length
      /** It will be stored decision SLAE */ 
      def xs = {val res = new Array[Double](n); Array.copy(b, 0, res, 0, n); res}
      // Philipp Ludwig von Seidel
      def xsSeidel = {solve(xs,xs); xs}
      def xsJacobi = {
        val xs2 = new Array[Double](n); solve(xs,xs2);
        xs
      }
      var METHOD_ID = 1
      /** error of solution */
      var ERROR = 1e-5
      
      /** @return solution of SLAE Ax=b (with ERROR) */
      def solve(xs1: Array[Double], xs2: Array[Double]): Unit = {
        /** norm(B) (see [1: p. 21, theor. 2]) */ 
        var q: Double = 0.0
        /** 
         * Необходимое количество итераций можно вычислить заранее!!! 
         * [1: 3.1, th. 2] -- {q = normM(B); dX = normV(x1-x0); e = ERROR;
         * eq: {(q^k)/(1-q)<=e*dX <=> q^k <= e*dX*(1-q), }} 
         **/
        var iter: Long = 0
        /**
         * see [1: p. 21, theor. 2]
         * return (norm(B)<1)
         */
        def isConverge : Boolean = {
          val diag = new Array[Double](n); for (i <- 0 until n) diag(i) = 1.0/A(i)(i);
          val B: Array[Array[Double]] = 
              if (METHOD_ID == 1)     //  Jacobi
                {
                  val res = new Array[Array[Double]](n);
                  for (i <- 0 until n) res(i) = Vector.BasicOperation.multi(A(i),-diag(i))
                  for (i <- 0 until n) res(i)(i) +=1
                  res
                } 
              else if (METHOD_ID == 2) // Seidel 
                {
                  val mat = Array.ofDim[Double](n,n);  
                  for (i <- 0 until n; j <- 0 to i) mat(i)(j) = -A(i)(j);
                  
                  val R = Array.ofDim[Double](n,n); 
                  for (i <- 0 until n; j <- i+1 until n) R(i)(j) = A(i)(j);
                  
                  val res = MatrixFunction.Multiply.multi(MatrixFunction.Inverse.inverseTriangleDown(mat),R)
                  res
                }
              else null
          q = MatrixFunction.Info.normMC(B);
          for (i <- 0 until n) diag(i)*=b(i)
          val xs2 = plus(MatrixFunction.Multiply.multi(B, xs), diag)
          val dX = Vector.Norm.norm1(minus(xs2, xs))
          iter = {
            val m1dX = 1.0/dX;
            var k: Long = 0; var mult = 1.0; 
            if (q<1.0) while (mult>ERROR*(1.0-q)*m1dX) {mult*=q; k+=1}
            else if (mult<ERROR*(1.0-q)*m1dX) {k = 1}
            else k = -1;
            k
          }
          Array.copy(xs2, 0, xs, 0, xs.length)
          val norm = MatrixFunction.Info.normMC(B)
          (norm<1.0)        
        }
        def isStability = {
          var flag  = true; 
          for (i <- 0 until n) flag&&=(!Number.isZero(A(i)(i))); 
          flag
        }
        val m1diagA = new Array[Double](n)
        for (i <- 0 until n) m1diagA(i) = 1.0/A(i)(i);
        
        val res = new Array[Double](n)
        // [1: 3.2]
        if (!isStability) throw new java.lang.IllegalArgumentException("На диагонали матрицы А есть нуль, решение невозможно")
        if (!isConverge) throw new java.lang.IllegalArgumentException("Решение невозможно -- норма матрицы B не меньше 1.0")
        if (iter == -1) throw new java.lang.IllegalArgumentException("Для решения необходимо слишком много итераций")
        while (iter>0){
          for (i <- 0 until n) 
            xs2(i) = -(Vector.BasicOperation.scalar(A(i),xs1) - A(i)(i)*xs1(i) - b(i))*m1diagA(i)
          iter-=1
        }
      }
}