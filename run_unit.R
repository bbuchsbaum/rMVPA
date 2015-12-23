library(RUnit)


testsuite <- defineTestSuite("rMVPA", dirs = "tests",
                             testFileRegexp = "^runit.+\\.R",
                             testFuncRegexp = "^test.+",
                             rngKind = "Marsaglia-Multicarry",
                             rngNormalKind = "Kinderman-Ramage")

testResult <- runTestSuite(testsuite)
printTextProtocol(testResult)
