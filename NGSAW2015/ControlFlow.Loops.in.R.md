## R语言中的控制流和循环

### 控制流

**标准的if else**

```r
p.test <- function(p) {
    if (p <= 0.05)
        print("yeah!!!!") else if (p >= 0.9)
        print("high!!!!") else print("somewhere in the middle")
}

p.test(0.5)
## [1] "somewhere in the middle"
```

**ifelse(test, yes, no)**

```r
p.test.2 <- function(p) {
    ifelse(p <= 0.05, print("yippee"), print("bummer, man"))
}

x <- runif(10, 0, 1)
x
##  [1] 0.27332 0.14155 0.89000 0.07041 0.79419 0.25013 0.02324 0.86766
##  [9] 0.41114 0.56165

p.test(x)

## Warning: the condition has length > 1 and only the first element will be used
## Warning: the condition has length > 1 and only the first element will be used

## [1] "somewhere in the middle"

p.test.2(x)

## [1] "yippee"
## [1] "bummer, man"

##  [1] "bummer, man" "bummer, man" "bummer, man" "bummer, man" "bummer, man"
##  [6] "bummer, man" "yippee"      "bummer, man" "bummer, man" "bummer, man"
```

### 其他的控制流向量化的方法

```r
p.1000 <- runif(n = 1000, min = 0, max = 1)
p.ifelse <- ifelse(p.1000 < 0.05, 1, 0)  # If it is less than 0.05, then you get a 1, otherwise 0.

sum(p.ifelse)/length(p.1000)
## [1] 0.059

#Same number, faster and simpler computation.
length(p.1000[p.1000 < 0.05])/length(p.1000)
## [1] 0.059
```

### 简单的loop

**while() {}**

```r
i <- 1
while (i <= 10) {
    print(i)
    i <- i + 0.5
}

## [1] 1
## [1] 1.5
## [1] 2
## [1] 2.5
## [1] 3
## [1] 3.5
## [1] 4
## [1] 4.5
## [1] 5
## [1] 5.5
## [1] 6
## [1] 6.5
## [1] 7
## [1] 7.5
## [1] 8
## [1] 8.5
## [1] 9
## [1] 9.5
## [1] 10
```

**for() {}**

```r

for (i in 1:10) {
    print(i)
}

## [1] 1
## [1] 2
## [1] 3
## [1] 4
## [1] 5
## [1] 6
## [1] 7
## [1] 8
## [1] 9
## [1] 10

for (i in seq(from = 1, to = 5, by = 0.5)) {
    print(i)
}

## [1] 1
## [1] 1.5
## [1] 2
## [1] 2.5
## [1] 3
## [1] 3.5
## [1] 4
## [1] 4.5
## [1] 5
```

另外，R中字符串无法进行如下操作：

```r
for (letter in “word”) {
    print(letter)
}

## [1] "word"
```

替代方法是：

```r

strsplit(“word”, split = “”)

## [[1]]
## [1] "w" "o" "r" "d"

for (letter in strsplit("word", split = "")) {
    print(letter)
}

## [1] "w" "o" "r" "d"
```

### 更多需要避免的循环

```r

for (i in 1:100) {
		print (rnorm(n=1, mean = 0, sd = 1))
}

## [1] -0.1837
## [1] -0.9313
## [1] 1.648
## [1] -0.6964
## [1] 0.2112
## [1] 0.3441
## [1] 1.036
## .........
## [1] -0.8865
## [1] 0.981
## [1] 0.5333
```

我们意在每个循环生成1个随机数，其实是不需要的，完全可以这样写：

```r
rnorm(n = 100, mean = 0, sd = 1)

##   [1] -0.08683 -1.55262 -1.16909  0.30451 -1.14555  0.76682  0.12643
##   [8] -0.61174 -0.29103 -0.10707 -0.03397 -0.05926  0.27294  1.32693
##  [15] -0.53284  1.83234  0.43959 -0.88991  0.25383  0.96709 -0.23210
##  [22] -1.00190 -1.32289  1.80030  1.15272 -1.82907  0.75989  1.35966
##  [29]  0.53943  0.01429 -0.58707 -0.11886 -0.70367 -2.38988  0.08033
##  [36] -0.22795 -0.62166 -0.19832 -1.95990 -0.85127  0.94236  0.37771
##  [43]  0.32617 -0.08393 -0.54506 -2.58781 -0.58433  0.20985 -0.41613
##  [50]  0.60527  0.51713  1.57950 -0.61079 -0.28564 -0.16444  0.55007
##  [57]  0.57258  0.58513 -0.86728 -0.81185 -0.29333 -1.23935  0.46169
##  [64] -1.53586 -0.32583  0.17629 -0.85579  1.04989  1.22120  1.53359
##  [71] -2.37276  1.44393  1.47506  0.40110 -0.10157  0.35485 -0.72068
##  [78] -1.27910  0.63152 -0.65216  1.60160  0.27109  0.50904 -1.00531
##  [85]  0.76743 -0.78954 -0.01159  1.06944  1.15661 -0.91031  1.54919
##  [92] -0.84334  2.19994  0.26716  0.02081  0.53577  0.07840 -0.79387
##  [99] -1.18941  1.24745
```

还有，当我们想要初始化一个向量以存储所有的数值时，会这样干：

```r
n <- 1e+5
head(x)

## [1] NA NA NA NA NA NA

for (i in 1:n) {
    x[i] <- rnorm(n = 1, mean = 0, sd = 1)
}
head(x)

## [1]  0.2848 -0.5432  1.1391 -1.0901  0.8515  0.5490
```

然而，这是十分低效的做法，很耗时间：

```r
system.time(

for (i in 1:n){
    x[i] <- rnorm(n=1, mean=0, sd=1)})

##    user  system elapsed
##   0.562   0.023   0.584

#We can also use the replicate function to do the same thing. Easier syntax to write.

system.time(z <- replicate(n, rnorm(n = 1, mean = 0, sd = 1)))

##    user  system elapsed
##   0.561   0.035   0.841
```

很显然，第二种做法也很低效。其实，R就是向量化的，可以直接这样做：

```r
x <- rnorm(n, 0, 1)
system.time(y <- rnorm(n, 0, 1))

##    user  system elapsed
##   0.010   0.000   0.011
```

因此，请记住，**在R中循环总是比函数家族的应用要慢（大数据除外），而后者要比向量化计算慢**。
