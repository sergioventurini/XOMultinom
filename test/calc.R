library(XOMultinom)

t <- 1
n <- 7
m <- 10
# debugonce(range_probability)
print(range_probability_C(t, n, m), digits = 15)

###

t <- 2
n <- 70
m <- 30
print(smallest_order_value_C(t, n, m), digits = 15)

###

t <- 6
n <- 60
m <- 30
J <- 3
# debugonce(highest_order_statistics)
print(highest_order_statistics_C(t, n, m, J), digits = 15)

###

t <- 2
n <- 6
m <- 3
print(max_order_statistic_C(t, n, m), digits = 15)

t <- 2
n <- 6
m <- 7
print(max_order_statistic_C(t, n, m), digits = 15)

t <- 32
n <- 60
m <- 3
print(max_order_statistic_C(t, n, m), digits = 15)

###

res <- numeric(15)

n <- 10
m <- 5

for (t in 2:16) {
  res[t - 1] <- max_order_statistic_C(t, n, m)
}
res

###

res <- numeric(19)

J <- 2

for (t in 5:23) {
  res[t - 4] <- highest_order_statistics_C(t, n, m, J)
}
res

###

res <- numeric(27)

n <- 15
m <- 15
J <- 3

for (t in 8:34) {
  res[t - 7] <- highest_order_statistics_C(t, n, m, J)
}
res

highest_order_statistics_C(8, 15, 15, 3)

###

res <- numeric(19)

for (t in 5:23) {
    res[t - 4] <- smallest_order_value_C(t, n, m)
}
res
