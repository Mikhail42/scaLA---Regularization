package LAfunction.Number;

import scala.math.abs
import scala.math.cos
/**
 * @author Ionkin Mikhail
 */
object DoubleFunction {
  import math._
  
  def sec(x: Double) = 1.0/cos(x)
  def two(x: Double) = x+x
  def four(x: Double) = two(x+x)
  
  def div2(x: Double) = x*0.5
  private val oneDiv3 = 1.0/3.0
  def div3(x: Double) = x*oneDiv3
  
  private var doubleError = 1e-10
  def setERROR(newERROR: Double) { doubleError = newERROR}
  
  def compareDouble(a: Double, b: Double) = 
    if (abs(a-b)<doubleError) 0 else if (a>b) 1 else -1
  
  def specSign(x: Double) = if (compareDouble(x,0.0)>=0) 1 else -1
  def isZero(x: Double) = abs(x)<doubleError
  /** x mod 2 == 0
   *  @param x is positive number
   */
  def isDiv2(x: Int) = (x&1)==0
}