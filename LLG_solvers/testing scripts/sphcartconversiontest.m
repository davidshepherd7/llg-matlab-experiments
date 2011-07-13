% Check that conversion is consistent for N random cartesian vectors
N = 100;
X_cart = 2*(rand(N,3) - 0.5);

for i=1:N
    X_sph(i,:) = carttosph(X_cart(i,:));
    
    X_cart2(i,:) = sphtocart(X_sph(i,:));
    
    result(i,:) = X_cart(i,:) - X_cart2(i,:);
end

max(result)