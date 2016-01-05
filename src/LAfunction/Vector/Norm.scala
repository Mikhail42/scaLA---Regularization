package LAfunction.Vector
import math._
import LAfunction.Vector._

  object Norm{
    /** max(abs(a(0:n-1))) */
    def normCubic(a: Array[Double], n: Int) = {
      var maxA = abs(a(0)); for (i <- 1 until n) maxA = max(maxA,abs(a(i))); maxA
    }
    /** normCubic(a,a.length) */
    def normCubic(a: Array[Double]): Double = normCubic(a,a.length);
    /** sum(abs(a(0:n-1))) */
    def norm1(a: Array[Double], n: Int): Double ={var sum = 0.0; for(i <- 0 until n) sum+=abs(a(i)); sum}
    /** norm1(a,a.length) */
    def norm1(a: Array[Double]): Double = norm1(a,a.length)
    /** sum(abs(A(0:n-1)(j))) */
    def norm1(A: Array[Array[Double]],j: Int, n: Int): Double ={var sum = 0.0; for(i <- 0 until n) sum+=abs(A(i)(j)); sum}
    /** norm1(A,j,A.length) */
    def norm1(A: Array[Array[Double]],j: Int): Double = norm1(A,j,A.length)
    /** sqrt(multi(a,a,0,n)) */
    def norm2(a: Array[Double], n: Int) = sqrt(LAfunction.Vector.BasicOperation.scalar(a,a,0,n))
    /** norm2(a,a.length) */
    def norm2(a: Array[Double]): Double = norm2(a,a.length)
  }