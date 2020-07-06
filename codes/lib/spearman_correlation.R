# spearman_rho_ci - 
# Input: 
#   - var1 
#   - var2 
#   - df 
# e.g. a dataframe (df) with two columns (var1 and var2)
# Returns Spearman rho and confidence interval 
spearman_rho_ci = function(var1, var2, df){
  cor = cor.test(var1, var2, 
                 method = 'spearman')
  R = cor$estimate
  n = nrow(df)
  z = 2
  a = 1.96
  L = ((1 + R- ((1 - R)* exp(2*z*(a/2)/sqrt(n - 3)))))/((1 + R + ((1 - R)* exp(2*z*(a/2)/sqrt(n - 3)))))
  U = ((1 + R- ((1 - R)* exp(-2*z*(a/2)/sqrt(n - 3)))))/((1 + R + ((1 - R)* exp(-2*z*(a/2)/sqrt(n - 3)))))
  return(c(cor$estimate, cor$p.value, L, U))
} 
# End spearman_rho_ci

