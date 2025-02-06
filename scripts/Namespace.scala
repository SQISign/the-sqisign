import io.StdIn.readLine

// 1. #define DISABLE_NAMESPACING in sqisign_namespace.h
// 2. build (cmake and make)
// 3. find . -name '*.a' -exec nm {} \; | grep '.c.o:\|T ' | scala ../scripts/Namespace.scala > sqisign_namespace.h
// 4. cp sqisign_namespace.h $SQISIGN_DIR/include

object Namespace extends App {

  val PREAMBLE = """
#ifndef SQISIGN_NAMESPACE_H
#define SQISIGN_NAMESPACE_H

//#define DISABLE_NAMESPACING

#if defined(_WIN32)
#define SQISIGN_API __declspec(dllexport)
#else
#define SQISIGN_API __attribute__((visibility("default")))
#endif

#define PARAM_JOIN3_(a, b, c) sqisign_##a##_##b##_##c
#define PARAM_JOIN3(a, b, c) PARAM_JOIN3_(a, b, c)
#define PARAM_NAME3(end, s) PARAM_JOIN3(SQISIGN_VARIANT, end, s)

#define PARAM_JOIN2_(a, b) sqisign_##a##_##b
#define PARAM_JOIN2(a, b) PARAM_JOIN2_(a, b)
#define PARAM_NAME2(end, s) PARAM_JOIN2(end, s)

#ifndef DISABLE_NAMESPACING
#define SQISIGN_NAMESPACE_GENERIC(s) PARAM_NAME2(gen, s)
#else
#define SQISIGN_NAMESPACE_GENERIC(s) s
#endif

#if defined(SQISIGN_VARIANT) && !defined(DISABLE_NAMESPACING)
#if defined(SQISIGN_BUILD_TYPE_REF)
#define SQISIGN_NAMESPACE(s) PARAM_NAME3(ref, s)
#elif defined(SQISIGN_BUILD_TYPE_OPT)
#define SQISIGN_NAMESPACE(s) PARAM_NAME3(opt, s)
#elif defined(SQISIGN_BUILD_TYPE_BROADWELL)
#define SQISIGN_NAMESPACE(s) PARAM_NAME3(broadwell, s)
#elif defined(SQISIGN_BUILD_TYPE_ARM64CRYPTO)
#define SQISIGN_NAMESPACE(s) PARAM_NAME3(arm64crypto, s)
#else
#error "Build type not known"
#endif

#else
#define SQISIGN_NAMESPACE(s) s
#endif
"""

  val EPILOGUE = """
#endif
"""

  val x = Iterator
    .continually(readLine)
    .takeWhile(_ != null).toList

  var scfile = ""
  val allFuns: List[(String, String)] = x.flatMap {
    case i if i.contains(".c.o:") =>
      scfile = i
      None
    case i =>
      i.split(" ").last match {
        case j if j.startsWith("_") =>
          Some((j.substring(1), scfile))
        case j =>
          Some((j, scfile))
      }
    
  }. // removing duplicates..
  groupBy(_._1).mapValues(k => k.distinct.toList.sortBy(_._2).reduceLeft((i,j) => ((i._1, s"${i._2}, ${j._2}")))).values.toList

  val maxFunLen = allFuns.map(i => i._1.length).max

  val filterFiles = List(
    "fips202.c",
    "tools.c",
    "randombytes_system.c",
    "randombytes_ctrdrbg.c",
    "randombytes_ctrdrbg_aesni.c",
    "foo.c",
    "aes_c.c",
    "aes_ni.c",
    "ctr_drbg.c"    
  )

  val genericFiles = List(
    // quaternion module
    "intbig.c",
    "algebra.c",
    "ideal.c",
    "dim4.c",
    "dim2.c",
    "integers.c",
    "lattice.c",
    "lat_ball.c",
    "finit.c",
    "printer.c",
    "rationals.c",
    "l2.c",
    "lll_verification.c",
    "lll_applications.c",
    "rationals.c",
    "normeq.c",
    "ibz_division.c",
    "hnf_internal.c",
    "hnf.c",
    "random_input_generation.c",
    "mem.c",
    // mp module
    "mp.c"
  ).map(i => s"$i.o:")

  val groupedByFile = 
    allFuns.
    groupBy(_._2).
    map(i => (i._1, i._2.distinct.sorted)).
    filter(i => filterFiles.forall(j => !i._1.contains(j))).toList.sortBy(_._1)

  println(PREAMBLE)



  groupedByFile.foreach(i => {
    println(s"// Namespacing symbols exported from ${i._1.replaceAll("\\.o:", "")}:")
    i._2.foreach(j => 
      println(s"#undef ${j._1}")
    )
    println
    i._2.foreach(j => {
      val padded = j._1.padTo(maxFunLen, " ").mkString
      if (genericFiles.contains(j._2)) {
        println(s"#define $padded SQISIGN_NAMESPACE_GENERIC(${j._1})")
      } else {
        println(s"#define $padded SQISIGN_NAMESPACE(${j._1})")
      }
    }
    )
    println
  })

  println(EPILOGUE)
}