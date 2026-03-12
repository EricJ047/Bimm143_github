# BIMM143 Class06:Function
Yuxuan Jiang PID:A17324184

- [Background](#background)
- [Our first function](#our-first-function)
- [A second function](#a-second-function)
- [A protein generating function](#a-protein-generating-function)

## Background

All functions in R have at least three things:

- A **name** that we use to call the function.
- One or more input **arguments**.
- The **body** the lines of R code that do the work.

## Our first function

Let’s write a silly little function called `add()` to add some numbers
(the input arguments).

``` r
add <- 
function(x,y) 
{
x+y
}
```

Now we can use this function

``` r
add(100,1)
```

    [1] 101

``` r
add(c(100,1,100),1)
```

    [1] 101   2 101

``` r
add(x=c(100,1,100),y=1)
```

    [1] 101   2 101

> Q. What if I give a multiple element vector to `x` and `y`?

``` r
add(x=c(100,1),y=c(100,1))
```

    [1] 200   2

> Q. What if I give three inputs to the function?

``` r
#add(x=c(100,1),y=1,z=1) error
```

> Q.What if I give only one input to the add function?

``` r
addnew <- 
function(x,y=1) 
{
x+y
}
```

``` r
addnew(x=100)
```

    [1] 101

``` r
addnew(c(100,1),100) 
```

    [1] 200 101

If we write our function with input arguments having no default value
then the user will be required to set them when they use the function.
We can give our input arguments “default” values by setting them equal
to some sensible value (e.g. y=1).

## A second function

Let’s try something more interesting: Make a sequence generating tool.

The `sample` function can be a useful starting point here:

``` r
sample(1:10,size=4)
```

    [1]  3  2  9 10

> Q.Generate 9 frandom numbers taken form the input vector x=1:10?

``` r
sample(1:10,size=9)
```

    [1]  3  8  1  9  7  2  5 10  4

> Q.Generate 12 frandom numbers taken form the input vector x=1:10?

``` r
sample(1:10,size=12, replace=T)
```

     [1] 1 1 7 8 3 4 8 2 7 9 5 5

> Q.Write code for the `sample()` function that generates nuceotide
> sequences of length 6?

``` r
sample(c('A','T','C','G'),size=6,replace=T)
```

    [1] "G" "T" "T" "T" "T" "G"

> Q.Write a first function `generate_dna()` that returns a user
> specified length DNA sequences:

``` r
generate_dna <- 
function(x=6) 
{
  sample(c('A','T','C','G'),size=x,replace=T)           
}
generate_dna(10)
```

     [1] "C" "T" "T" "G" "T" "T" "A" "C" "G" "T"

> **Key-Points** Every function in R looks fundamentally the same in
> terms of its structure. Basically 3 things:name, input, and body.

    function_name <- 
    function(input)
    {
    body
    }

> Functions can have multiple inputs. These can be **required**
> arguments or **optional** arguments, with optional arguments having a
> set of defult value.

> Q.Modify and improve our `generate_dna()` function to return it’s
> generated sequence in a more standard format like “AGTAGTA” rather
> than the vector “A”,“G”,“T”,“A”.

``` r
generate_dna <- 
function(x=6,fasta=T) 
{
 ans <- sample(c('A','T','C','G'),size=x,replace=T)
 if(fasta)
 {
    cat("Single-element vector output")
    ans <- paste(ans,collapse = "")
 }else
  {
    cat("Multi-element vector output")
  }
 return(ans)
}
generate_dna(10,T)
```

    Single-element vector output

    [1] "GGGATGTAAA"

``` r
generate_dna(10,F)
```

    Multi-element vector output

     [1] "C" "C" "C" "T" "G" "T" "G" "T" "G" "C"

The `paste()` function - it’s job is to join up or stick together
(aka.paste) input strings together.

``` r
paste("alice","loves R",sep='****')
```

    [1] "alice****loves R"

Flow control means where the R brain goes in your code.

``` r
good_mood <- T
if(good_mood)
{
  cat("Great!")
}else
 {
   cat("Bummer!")
 }
```

    Great!

## A protein generating function

> Q.Write a function that generates a user specified length protein
> sequence.

``` r
# The aa to sample from
aabank <-c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V")
gen_protein <- 
function(x)
{
# Draw n=x aa to make our sequence
  protein <- sample(aabank,size=x,replace=T)
  protein <- paste(protein,collapse = "")
  return (protein)
}
```

> Q.Use that function to generate random protein sequences between
> length 6 and 12.

``` r
for(i in 6:12)
{
  #FASTA ID line
  cat(">",i,sep="","\n")
  #Sequecen line
  cat(gen_protein(i),"\n")
}
```

    >6
    GRHHAI 
    >7
    WCLVWMP 
    >8
    SPMMACYP 
    >9
    LNRVQFKIS 
    >10
    EVGMVLGNMK 
    >11
    PPMPEHGRAKP 
    >12
    LRMQEILSCESR 

> Q.Are any of your sequences unique i.e.not found anywhere in nature?

Generated sequences with length\>=8 are unique based on the protein
BLAST result.
