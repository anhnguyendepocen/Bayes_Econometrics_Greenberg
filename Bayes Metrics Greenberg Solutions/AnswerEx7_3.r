library(MCMCpack)
"CoalDisast" <-
structure(list(Year = as.integer(c(1851, 1852, 1853, 1854, 1855, 
1856, 1857, 1858, 1859, 1860, 1861, 1862, 1863, 1864, 1865, 1866, 
1867, 1868, 1869, 1870, 1871, 1872, 1873, 1874, 1875, 1876, 1877, 
1878, 1879, 1880, 1881, 1882, 1883, 1884, 1885, 1886, 1887, 1888, 
1889, 1890, 1891, 1892, 1893, 1894, 1895, 1896, 1897, 1898, 1899, 
1900, 1901, 1902, 1903, 1904, 1905, 1906, 1907, 1908, 1909, 1910, 
1911, 1912, 1913, 1914, 1915, 1916, 1917, 1918, 1919, 1920, 1921, 
1922, 1923, 1924, 1925, 1926, 1927, 1928, 1929, 1930, 1931, 1932, 
1933, 1934, 1935, 1936, 1937, 1938, 1939, 1940, 1941, 1942, 1943, 
1944, 1945, 1946, 1947, 1948, 1949, 1950, 1951, 1952, 1953, 1954, 
1955, 1956, 1957, 1958, 1959, 1960, 1961, 1962)), Count = c(4, 
5, 4, 1, 0, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6, 3, 3, 5, 4, 5, 3, 
1, 4, 4, 1, 5, 5, 3, 4, 2, 5, 2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 
1, 1, 3, 0, 0, 1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1, 
0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2, 3, 3, 1, 1, 2, 
1, 1, 1, 1, 2, 4, 2, 0, 0, 0, 1, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
1, 0, 0, 1, 0, 1)), .Names = c("Year", "Count"), row.names = c("1", 
"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
"14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", 
"25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", 
"36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", 
"47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", 
"58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", 
"69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", 
"80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", 
"91", "92", "93", "94", "95", "96", "97", "98", "99", "100", 
"101", "102", "103", "104", "105", "106", "107", "108", "109", 
"110", "111", "112"), class = "data.frame")


windows(record = T)
y <- as.matrix(CoalDisast$Count, nrow = 112)

poissonModel1.out <- MCMCpoissonChange(y ~ 1, m = 1, c0 = 8, d0 = 1, mcmc = 10000,
     marginal.likelihood = c("Chib95"))
print("Summary:  One change point.")	 
summary(poissonModel1.out)	 
plot(poissonModel1.out)

poissonModel0.out <- MCMCpoisson(y ~ 1, c0 = 8, d0 = 1, mcmc = 10000, B0 = 1,
     marginal.likelihood = c("Laplace"))
	 	 
	 
print("Summary:  No change point.")	 
summary(poissonModel0.out)

#  Convert MLs from base e to base 10.
logMLOneChange <- format(log10(exp(attr(poissonModel1.out, "logmarglike"))), digits = 5)	 
paste("log base 10 marg like, One change = ", logMLOneChange)

logMLNoChange <- format(log10(exp(attr(poissonModel0.out, "logmarglike"))), digits = 5)	 
paste("log base 10 marg like, No change = ", logMLNoChange)

op <- par()
par(mfrow=c(attr(poissonModel1.out, "m") + 1, 1), mai=c(0.4, 0.6, 0.3, 0.05))
plotState(poissonModel1.out, legend.control = c(1, 0.6))
plotChangepoint(poissonModel1.out, verbose = TRUE, ylab="Density", start=1, overlay=TRUE)
par(op)	