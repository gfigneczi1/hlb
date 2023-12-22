function calcL2NormOfPMatrices(parameters)
fn = fieldnames(parameters);
for i=1:length(fn)
    for j=1:length(fn)
        A = parameters.(fn{i});
        A = reshape(A,3,7);
        B = parameters.(fn{j});
        B = reshape(B,3,7);
        C = A-B;
        n(i,j)=sqrt( max( eig(C.'*C) ) );
    end
end
n = n/max(max(n));
m = (n<0.2).*(n>0);


for i=1:length(fn)
    for j=1:length(fn)
        A = parameters.(fn{i});
        A = reshape(A,3,7);
        B = parameters.(fn{j});
        B = reshape(B,3,7);
        C = A(:,1:3)-B(:,1:3);
        n(i,j)=sqrt( max( eig(C.'*C) ) );
    end
end

for i=1:length(fn)
    for j=1:length(fn)
        A = parameters.(fn{i});
        A = reshape(A,3,7);
        B = parameters.(fn{j});
        B = reshape(B,3,7);
        C = A(:,4:6)-B(:,4:6);
        n(i,j)=sqrt( max( eig(C.'*C) ) );
    end
end

for i=1:length(fn)
    for j=1:length(fn)
        A = parameters.(fn{i});
        A = reshape(A,3,7);
        B = parameters.(fn{j});
        B = reshape(B,3,7);
        C = A(:,7)-B(:,7);
        n(i,j)=sqrt( max( eig(C.'*C) ) );
    end
end

end

