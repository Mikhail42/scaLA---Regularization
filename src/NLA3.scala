import math._
import LAfunction._
import LAfunction.Number.DoubleFunction._

object NLA3 {
  import Data._
  import BasicData._
  import Cond._
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
  
  /** \int_{a}^{b} {K(s,t)\phi(s)} \dd[s] = f(t) */
  def main(args: Array[String]): Unit = {
    Data.BasicData.a
    println("A = ")
    println(LAfunction.MatrixFunction.IO.toString(A))
    println("cond_2(A)     = "+condA._1
      +(if (condA._2==null) "" else ". cond_2(A) не точно"))
    println("cond_2(A^T A) = "+condATA._1
      +(if (condATA._2==null) "" else ". cond_2(A^T A) не точно"))
    println
    for (i <- 0 until Alpha.nBestSqrAlphas)
      println("cond_2(A^T A + diag[{"+bestSqrAlphas(i)+"}]) = "+vectCond(i))
    println
    println("Числа обусловленности матрицы")
    println("    "+"(\u03B2E      A   )")
    println("A = "+"(            )")
    println("    "+"(A^T (\u03B1^2/\u03B2)E)")
    println("при \u03B1^2 = "+LAfunction.Vector.Print.toString(Data.Alpha.bestSqrAlphas, 7)+", ")
    println("при \u03B2 = {1, sqrt(\u03B1), sqrt(\u03B1^2+0.5*\u03C3_m^2(A))}")
    println(LAfunction.MatrixFunction.IO.toString(Data.Cond.matCond))
    println    
    
    val xs = Data.BasicData.xs
    
    val frames = new Array[(draw.MyFrame, String)](4)
    
    val luFrame = new draw.MyFrame; frames(0) = (luFrame, "LU")
    val qrFrame = new draw.MyFrame; frames(1) = (qrFrame, "QR")
    val lupFrame = new draw.MyFrame; frames(3) = (lupFrame, "LUP")
    val GSFrame = new draw.MyFrame; frames(2) = (GSFrame, "Gauss-Seidel")
    1
    for (frame <- frames) {
      try{
        for (i <- 0 until Data.Alpha.nBestSqrAlphas){
          val solut =  frame._2 match {
            case "LU" => Data.Solution.LU.solveATAsWithAlpha(i)
            case "LUP" => Data.Solution.LUP.solveATAsWithAlpha(i)
            case "Gauss-Seidel" => Data.Solution.Jacobi.solveATAsWithAlpha(i)
            case _ => Data.Solution.QR.solveATAsWithAlpha(i)
          } 
          frame._1.addGraph(xs, solut, "\u03B1^2="+Data.Alpha.bestSqrAlphas(i))
        }
        for (i <- 0 until Data.Beta.nBetas){
          val solut =  frame._2 match {
            case "LU" => Data.Solution.LU.solveATAsWithAlphaAndBeta(i)
            case "LUP" => Data.Solution.LUP.solveATAsWithAlphaAndBeta(i)
            case "Gauss-Seidel" => Data.Solution.Jacobi.solveATAsWithAlphaAndBeta(i)
            case _ => Data.Solution.QR.solveATAsWithAlphaAndBeta(i)
          }        
          frame._1.addGraph(xs, solut, "\u03B1="+Alpha.bestSqrAlphas(3) +", \u03B2="+Beta.namesBeta(i))
        }
        frame._1.building(frame._2+": "+Data.BasicData.name)
      } catch {
        case e: Exception => {}
      }
    }
  }
}