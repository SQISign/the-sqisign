/**
 * Code to analyze instrumented code from the SQIsign IntBig module.
 *
 * Features:
 * - verifies arithmetic
 * - aggregate number of errors / ok per function
 * - aggregate minimum / maximum values per function
 *
 * Prerequisite: enable debug output in intbig.x: #define DEBUG_VERBOSE
 * Usage: ./<test app> | scala IntBigTest.scala
 * Usage with unit test: sqisign_test_intbig <reps> <bits> | scala IntBigTest.scala
 *
 * Run option -v: verbose full output
 */

object IntBigTest {

  // Test functions
  object IntBigTestFuns {
    def ibz_add(a: Array[BigInt]) = IntBigRes(a(0) == a(1) + a(2), a(1) + a(2), a)
    def ibz_sub(a: Array[BigInt]) = IntBigRes(a(0) == a(1) - a(2), a(1) - a(2), a)
    def ibz_mul(a: Array[BigInt]) = IntBigRes(a(0) == a(1) * a(2), a(1) * a(2), a)
    def ibz_div(a: Array[BigInt]) = IntBigRes(a(0) == a(2) / a(3) && a(1) == a(2) % a(3), a(2) / a(3), a)
    def ibz_pow_mod(a: Array[BigInt]) = IntBigRes(a(0) == a(1).modPow(a(2), a(3)),  a(1).modPow(a(2), a(3)), a)
    def ibz_cmp(a: Array[BigInt]) = IntBigRes(
      if (a(1) == a(2)) a(0) == 0 else if (a(1) < a(2)) a(0) < 0 else a(0) > 0,
      if (a(1) == a(2)) -1 else if (a(1) < a(2)) 1 else 0,
      a)
    def ibz_is_zero(a: Array[BigInt]) = IntBigRes(if (a(1) == 0) a(0) == 1 else a(0) == 0, if (a(1) == 0) 1 else 0, a)
    def ibz_is_one(a: Array[BigInt]) = IntBigRes(if (a(1) == 1) a(0) == 1 else a(0) == 0, if (a(1) == 1) 1 else 0, a)
    def ibz_probab_prime(a: Array[BigInt]) = IntBigRes(if (a(1).isProbablePrime(a(2).toInt)) a(0) > 0 else a(0) == 0, if (a(1).isProbablePrime(a(2).toInt)) 1 else 0, a)
    def ibz_gcd(a: Array[BigInt]) = IntBigRes(a(1).gcd(a(2)) == a(0), a(1).gcd(a(2)), a)
    def ibz_sqrt_mod_p(in: Array[BigInt]): IntBigRes = {
      val sqrt = in(0)
      val p = in(2)
      val a = if (in(1).mod(p) < 0) p + in(1).mod(p) else in(1).mod(p)
      val exp0 = sqrt.modPow(2, p)
      IntBigRes(exp0 == a || (p - exp0) == a, sqrt.modPow(2, p), in)
    }
    def ibz_sqrt_mod_2p(a: Array[BigInt]): IntBigRes = IntBigRes(a(0).modPow(2, 2 * a(2)) == a(1), a(0), a)
  }

  val funList = Map(
    //"ibz_add" -> ibz_add _,
    "ibz_sqrt_mod_p" -> IntBigTestFuns.ibz_sqrt_mod_p _,
    "ibz_sqrt_mod_2p" -> IntBigTestFuns.ibz_sqrt_mod_2p _,
    "ibz_add" -> IntBigTestFuns.ibz_add _,
    "ibz_sub" -> IntBigTestFuns.ibz_sub _,
    "ibz_mul" -> IntBigTestFuns.ibz_mul _,
    "ibz_div" -> IntBigTestFuns.ibz_div _,
    "ibz_pow_mod" -> IntBigTestFuns.ibz_pow_mod _,
    "ibz_cmp" -> IntBigTestFuns.ibz_cmp _,
    "ibz_is_zero" -> IntBigTestFuns.ibz_is_zero _,
    "ibz_is_one" -> IntBigTestFuns.ibz_is_one _,
    "ibz_probab_prime" -> IntBigTestFuns.ibz_probab_prime _,
    "ibz_gcd" -> IntBigTestFuns.ibz_gcd _
  )

  case class AggregateResults(funName: String, errors: Int, ok: Int, max: Option[Int], min: Option[Int]) {
    def errInc = AggregateResults(funName, errors + 1, ok, max, min)
    def okInc(operands: List[BigInt]) = {
      val operandsBitLen = operands.map(_.bitLength)
      val newMax = Some((max.getOrElse(0) :: operandsBitLen).max)
      val newMin = Some((min.getOrElse(Int.MaxValue) :: operandsBitLen).min)
      AggregateResults(funName, errors, ok + 1, newMax, newMin)
    }

    override def toString: String = s"$funName: $errors errors, $ok ok, max value: ${max.getOrElse(BigInt(0))} bits, min value: ${min.getOrElse(BigInt(0))} bits)"
  }

  def main(args: Array[String]) = {
    val v = args.length >= 1 && args(0) == "-v"
    var err = Map() ++ funList.map(i => (i._1 -> AggregateResults(i._1, 0, 0, None, None)))
    var cont = true
    while (cont) {
      val l = scala.io.StdIn.readLine()
      if (l == null) {
        val numerr =
        err.foreach(i => println(i._2))
        println(s"==========\n${err.values.map(_.errors).sum} errors found\n${err.values.map(_.ok).sum} checks ok")
        cont = false
      } else {
        err = check(l, v, err)
      }
    }
  }

  case class IntBigRes(verif: Boolean, expected: BigInt, got: Array[BigInt])

  def check(line: String, v: Boolean, agg: Map[String, AggregateResults]): Map[String, AggregateResults] = {
    line.split(",").toList match {
      case x :: xs if funList.contains(x) =>
        val funA = xs.map(i => BigInt(i, 16)).toArray
        funList(x)(funA) match {
          case IntBigRes(false, exp, got) =>
            println(s"function: $x\ngot:\n${got.map(_.toString(16)).mkString(",")}\nexpected:\n${exp.toString(16)}")
            agg.map(i => if (i._1 == x) i._1 -> i._2.errInc else i)
          case _ =>
            agg.map(i => if (i._1 == x) i._1 -> i._2.okInc(funA.toList) else i)
        }
      case _ =>
        if (v) println(line)
        agg
    }
  }

}
