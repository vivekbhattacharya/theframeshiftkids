### Replacing if-loops ###
Instead of
```
for i = 1:length(v)
    if v(i) == 2
       v(i) = 3;
    end
end
```
do
```
x = find(v == 2);
v(x) = 3;
```
Because the latter throws all the hard work to Matlab's compiled internal code, it'll be much faster. It works because logical tests return arrays. (Try chaining together more boolean operators.)

### Generating ranges in matrix format ###
How do you get 1:84 into a 28x3 matrix? `reshape(1:84, 3, 28)'`. Note that the apostrophe at the end transposes the matrix because `reshape` rejiggers your matrices in one very specific way: up to down first, then left to right). For more information, consult `doc reshape`.