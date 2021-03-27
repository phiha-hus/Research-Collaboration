function B = binom(I,J)

for i=0:I
    for j=0:J
        B(i+1,j+1)=faculty(i)/(faculty(i-j)*faculty(j));
    end
end