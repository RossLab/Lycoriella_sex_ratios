# October 2024
# Lycoriella ingenua sex ratio experiment

library(ggplot2)
library(dplyr)
library(tidyverse)
library(data.table)
library(lme4)
library(stats)
library(lmerTest)
library(sjPlot)
library(meta)
library(nnet)
library(broom)
library(emmeans)
library(car)

setwd('/Users/robertbaird/Library/Mobile Documents/com~apple~CloudDocs/Documents/analyses/fungus_gnats/Ling_SR/')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# LOADING IN DATA

# sex_ratios_all file - sex ratios for every single vial in this study (excluding 12 and 25 temperature treatments)
sex_ratios_all <- read.table('Ling_SR_data_all_vials.txt', header=T, stringsAsFactors=F, fill=T)
head(sex_ratios_all)
nrow(sex_ratios_all)
# sex_ratios_MD - sex ratios for every mother-daughter pair - has fewer rows than file above because the G0 daughter vials don't have a mother record to be paired with 
sex_ratios_MD <- read.table('Ling_mother_daughter_sex_ratios.txt', header=T, stringsAsFactors=F, fill=T)
head(sex_ratios_MD)
nrow(sex_ratios_MD)
# temperature experiment results
temp <- read.table('Ling_SR_data_temperatures.txt', header=T, stringsAsFactors=F, fill=T)
head(temp)
nrow(temp)

# Some of the main questions we want to answer:

# 1. Are primary sex ratios variable in this species, i.e. do they significantly deviate from 1:1? 
# 2. Is the primary sex ratio heritable: is there a correlation between the primary sex ratio that a mother produces and the primary sex ratio that her daughter produces? Do siblings produce more similar primary sex ratios than non-sibs? Do females within an isofemale line produce more similar primary sex ratios?
# 3. Is there a change in the primary sex ratios over generations? Do they become more skewed, or settle on 1:1?
# 4. Are any of the above changes due to differences in female vs male mortality?
# 5. Is there an effect of temperature on the sex ratio?

# Bear in mind that some brood sizes were very low (e.g. some were as low as 1), while the highest was 109. I don't want flies that had extremely low counts to skew the result, i.e. they should be to some extent weighted by the total brood size.

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

### 1. Testing if primary sex ratios deviate from 1:1 (Binomial Test)

# Approach: 1. Perform binomial tests for each brood, 2. Use a Fisher’s method, a meta-analysis technique to combine p-values. For each brood, conduct a binomial test, where the null hypothesis is that the probability of producing a male or female is 0.5.

# Add total brood size column
sex_ratios_all$brood_size <- sex_ratios_all$males+sex_ratios_all$females
sex_ratios_all <- sex_ratios_all[complete.cases(sex_ratios_all),]

# Check means and medians
#> median(sex_ratios_all$sex_ratio)
#[1] 0.5365854
#> mean(sex_ratios_all$sex_ratio)
#[1] 0.5218199

# Binomial test for each brood
sex_ratios_all <- sex_ratios_all %>%
  rowwise() %>%
  mutate(p_value = binom.test(males, brood_size, p = 0.5)$p.value)

# Output the results
print(sex_ratios_all)
# Visualize the distribution of p-values
hist(p_values, main = "Histogram of Individual Fisher's Test p-values", xlab = "p-values")

# Fisher's method combines p-values by transforming them into test statistics.
p_values <- sex_ratios_all$p_value  # vector of p-values from individual tests
# Compute the test statistic for Fisher's method using the negative log of p-values.
# Calculate Fisher's test statistic
fisher_stat <- -2 * sum(log(p_values))
# Calculate degrees of freedom (2 times the number of p-values)
df_fisher <- 2 * length(p_values)
# Get the combined p-value from the chi-squared distribution
fisher_p_value <- pchisq(fisher_stat, df = df_fisher, lower.tail = FALSE)
# Output the result
cat("Fisher's combined p-value:", fisher_p_value, "\n") # 3.322595e-136 

# Alternative meta-analysis approach: Stouffer’s method

# Calculate Z-scores for individual p-values
z_scores <- qnorm(1 - p_values / 2)
# Calculate weights based on the total offspring count
total_counts <- sex_ratios_all$brood_size
weights <- sqrt(total_counts)
# Combine Z-scores using Stouffer's method
combined_z <- sum(weights * z_scores) / sqrt(sum(weights^2))
# Convert combined Z-score back to a p-value
combined_p_value_stouffer <- 2 * (1 - pnorm(abs(combined_z)))
cat("Stouffer's combined p-value:", combined_p_value_stouffer, "\n") # 0

# Plotting - density plot

p <- ggplot(sex_ratios_all, aes(x=sex_ratio, fill=generation)) +
  geom_density(alpha=.25) +
  theme_classic() +
  xlab("Proportion of male offspring") +
  ylab("Density") +
  scale_fill_discrete(name="Generation", breaks=c('G0', 'F1', 'F2', 'F3', 'F4')) +
  theme(axis.line=element_line(linewidth=1),
        axis.ticks=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black", face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=14, face="bold"))

print(p)
#ggsave(file="/Users/robertbaird/documents/writing/Ling_SR/Figure_bits/SR_density.svg", plot=p, width=5.5, height=4)#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

### 2. Testing correlation of sex ratios between mothers and daughters (regression/correlation)

# Approach: To test the heritability of the primary sex ratio, calculate the sex ratio for each mother and her daughter and then use e.g. Pearson correlation or linear regression to determine if there is a significant relationship.

head(sex_ratios_MD)
sex_ratios_MD <- sex_ratios_MD[complete.cases(sex_ratios_MD),]
# add total brood size column
sex_ratios_MD$mother_brood_size <- sex_ratios_MD$mother_males+sex_ratios_MD$mother_females
sex_ratios_MD$daughter_brood_size <- sex_ratios_MD$daughter_males+sex_ratios_MD$daughter_females
# add ratio column
sex_ratios_MD$mother_ratio <- sex_ratios_MD$mother_males/sex_ratios_MD$mother_brood_size
sex_ratios_MD$daughter_ratio <- sex_ratios_MD$daughter_males/sex_ratios_MD$daughter_brood_size

# Create new columns where both mother_ratio and daughter_ratio are weighted by their respective brood sizes
sex_ratios_MD$weighted_mother_ratio <- sex_ratios_MD$mother_ratio * sqrt(sex_ratios_MD$mother_brood_size)
sex_ratios_MD$weighted_daughter_ratio <- sex_ratios_MD$daughter_ratio * sqrt(sex_ratios_MD$daughter_brood_size)

# Weighted linear regression using the new weighted columns
model <- lm(weighted_daughter_ratio ~ weighted_mother_ratio, data = sex_ratios_MD)

# Summary of the model
summary(model)

#Call:
#  lm(formula = weighted_daughter_ratio ~ weighted_mother_ratio, 
#     data = sex_ratios_MD)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-3.2911 -1.0000 -0.1328  0.8800  5.8095 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            1.49178    0.15001   9.945  < 2e-16 ***
#  weighted_mother_ratio  0.26594    0.04009   6.633 9.74e-11 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 1.45 on 438 degrees of freedom
#Multiple R-squared:  0.09127,	Adjusted R-squared:  0.0892
#F-statistic: 43.99 on 1 and 438 DF,  p-value: 9.737e-11


## can also see  if siblings (daughters sharing the same mother_vial_id) have more similar sex ratios than non-siblings.
# do this by :
# 1.	Calculate Absolute Differences: Calculate the absolute differences in sex ratios for both siblings and non-siblings.
# 2.	Weighted Differences: Adjust the differences using brood sizes to ensure that smaller broods do not disproportionately affect the analysis.
# 3.	Statistical Testing: Perform statistical tests to compare the similarity in sex ratios between the two groups.

# Calculate absolute differences for siblings
sibling_diffs <- sex_ratios_MD %>%
  inner_join(sex_ratios_MD, by = "mother_vial_id", suffix = c(".sibling", ".sibling_pair")) %>%
  filter(daughter_vial_id.sibling != daughter_vial_id.sibling_pair) %>%
  mutate(diff = abs(daughter_ratio.sibling - daughter_ratio.sibling_pair),
         weighted_diff = diff * sqrt(daughter_brood_size.sibling * daughter_brood_size.sibling_pair)) %>%
  select(daughter_vial_id.sibling, daughter_vial_id.sibling_pair, diff, weighted_diff)

# Calculate absolute differences for non-siblings
non_sibling_diffs <- sex_ratios_MD %>%
  inner_join(sex_ratios_MD, by = character(), suffix = c(".non_sibling", ".non_sibling_pair")) %>%
  filter(mother_vial_id.non_sibling != mother_vial_id.non_sibling_pair) %>%
  mutate(diff = abs(daughter_ratio.non_sibling - daughter_ratio.non_sibling_pair),
         weighted_diff = diff * sqrt(daughter_brood_size.non_sibling * daughter_brood_size.non_sibling_pair)) %>%
  select(daughter_vial_id.non_sibling, daughter_vial_id.non_sibling_pair, diff, weighted_diff)

# Display the first few rows of sibling differences
head(sibling_diffs)
# Display the first few rows of non-sibling differences
head(non_sibling_diffs)

# Calculate means and standard deviations for siblings
sibling_stats <- sibling_diffs %>%
  summarise(mean_diff = mean(weighted_diff),
            sd_diff = sd(weighted_diff),
            n = n())

# Calculate means and standard deviations for non-siblings
non_sibling_stats <- non_sibling_diffs %>%
  summarise(mean_diff = mean(weighted_diff),
            sd_diff = sd(weighted_diff),
            n = n())

# Print results
print(sibling_stats)
print(non_sibling_stats)

# Perform t-test on weighted differences
t_test_result <- t.test(sibling_diffs$weighted_diff, non_sibling_diffs$weighted_diff)

# Print t-test results
print(t_test_result)
# t = -9.4917, df = 1892.5, p-value < 2.2e-16

# Plotting mother-daughter correlation

# exclude brood sizes smaller than 10 so it looks less messy
sex_ratios_plot <- sex_ratios_MD[which(sex_ratios_MD$daughter_brood_size > 10),]

p <- ggplot(sex_ratios_plot, aes(mother_ratio, daughter_ratio))+
  geom_point(aes(size=daughter_brood_size), alpha=0.7)+
  geom_smooth(method = "lm", se=F, colour="black")+
  theme_classic()+
  xlim(0,1)+
  ylim(0,1)+
  geom_vline(xintercept=c(0.5), linetype="dotted", size=1, colour="black")+
  geom_hline(yintercept=c(0.5), linetype="dotted", size=1, colour="black")+
  xlab("Proportion of male offspring (mother)")+
  ylab("Proportion of male offspring (daughter)")+
  theme(axis.line=element_line(linewidth=1),
        axis.ticks=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black", face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=14, face="bold")) +
  labs(size="Daughter Brood Size")
print(p)
#ggsave(file="/Users/robertbaird/documents/writing/Ling_SR/Figure_bits/mother_daughter_SRs.svg", plot=p, width=7.5, height=5.2)#

# seeing if isofemale lines differ from one another in their variance

head(sex_ratios_MD)

# Levene's test (for homogeneity of variance)
levene_test <- leveneTest(daughter_ratio ~ founder_id, data = sex_ratios_MD)
print(levene_test)

# plot sex ratios within each isofemale line for supplementary figure
p <- ggplot(sex_ratios_MD, aes(mother_ratio, daughter_ratio)) +
  geom_point(aes(size = daughter_brood_size), alpha = 0.7) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.5))  +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.5))  +
  geom_vline(xintercept = 0.5, linetype = "dotted", size = 1, colour = "black") +
  geom_hline(yintercept = 0.5, linetype = "dotted", size = 1, colour = "black") +
  xlab("Proportion of male offspring (mother)") +
  ylab("Proportion of male offspring (daughter)") +
  theme(
    axis.ticks = element_line(linewidth = 1),        # Keep axis ticks for left and bottom
    axis.text = element_text(size = 14, colour = "black", face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    strip.text = element_blank(),    # Remove the facet labels (founder_id)
    strip.background = element_blank()  # Remove the strip background
  ) +
  labs(size = "Daughter Brood Size") +
  facet_wrap(~founder_id, ncol = 7, nrow = 4) +
  theme(
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.spacing = unit(1, "lines"),               # Adjust spacing between panels
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1)  # Keep panel borders
  )
print(p)

ggsave(file="/Users/robertbaird/documents/writing/Ling_SR/Figure_bits/mother_vs_daughter_by_line_revised.svg", plot=p, width=12.5, height=6)#


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

### 3. Testing changes in sex ratio across generations (Mixed-effects model)

# Approach: To assess whether primary sex ratios change over generations using a mixed-effects model, with sex ratio as the response variable and generation as a fixed effect. A mixed-effects model can account for the fact that broods are nested within families and handle the repeated measures across generations.

head(sex_ratios_all)

#sex_ratios_all <- within(sex_ratios_all, generation[generation == "G0"] <- 'F0') # for pairwise comparisons

# add sex ratio column
sex_ratios_all$sex_ratio <- sex_ratios_all$males/sex_ratios_all$brood_size

# Define bias categories based on sex ratio
sex_ratios_all$bias_category <- ifelse(sex_ratios_all$sex_ratio > 0.6, "Male-biased",
                           ifelse(sex_ratios_all$sex_ratio < 0.4, "Female-biased", "Unbiased"))

# Check the table
table(sex_ratios_all$generation, sex_ratios_all$bias_category)

# Fit the weighted multinomial logistic regression model
model_multinom_weighted <- multinom(bias_category ~ generation, data = sex_ratios_all, weights = sex_ratios_all$brood_size)

# Summary of the model
summary(model_multinom_weighted)

# pairwise comparisons between generations
# compute estimated marginal means for the bias_category by generation
emmeans_results <- emmeans(model_multinom_weighted, ~ generation | bias_category)
# pairwise comparisons between the generations
pairwise_comparisons <- pairs(emmeans_results)
summary(pairwise_comparisons)

# Plot the proportion of bias categories over generations
sex_ratios_all$generation <- factor(sex_ratios_all$generation, levels = c("G0", "F1", "F2", "F3", "F4", "F5"))
sex_ratios_all$bias_category <- factor(sex_ratios_all$bias_category, levels = c("Male-biased", "Unbiased", "Female-biased"))
p <- ggplot(sex_ratios_all, aes(x = factor(generation), fill = factor(bias_category), weight = brood_size)) +
  geom_bar(position = "fill", colour="black") +
  theme_classic() +
  labs(x = "Generation", y = "Weighted Proportion") +
  labs(fill = "Bias category") +
  scale_fill_manual(values = c("Male-biased" = "royalblue3", "Unbiased" = "gray", "Female-biased" = "orange")) +  # Color for male and female boxes
  theme(axis.line=element_line(linewidth=1),
        axis.ticks=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black", face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=14, face="bold"))

print(p)

#ggsave(file="/Users/robertbaird/documents/writing/Ling_SR/Figure_bits/change_over_generations.svg", plot=p, width=5.5, height=4)#

# Key Takeaways:
# The log-odds of broods being either Male-biased or Unbiased decrease significantly over generations F3 and F4 compared to being Female-biased.
# This suggests that over time, there are fewer male-biased and Unbiased broods, and more broods tend to become female-biased.

# but is this due to an increase in male mortality...? (next question)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

### 4. Testing for a decrease in clutch size over generations

#use a linear regression or a mixed-effects model to assess whether clutch size (total offspring) decreases over generations

head(sex_ratios_all)

#sex_ratios_all <- within(sex_ratios_all, generation[generation == "G0"] <- 'F0') # for pairwise comparisons

# Mixed-effects model with clutch size (total_brood_size) as the response and generation as fixed effect
clutch_model <- lmer(brood_size ~ generation + (1 | founder_id), data = sex_ratios_all)

# Summary of the model
summary(clutch_model)

# Use emmeans to get pairwise contrasts
generation_emmeans <- emmeans(clutch_model, ~ generation)

# Get p-values for the specific contrasts of interest
pairwise_contrasts <- contrast(generation_emmeans, method = "consec") # "consec" tests consecutive levels

# Display the contrasts and their p-values
summary(pairwise_contrasts)

#contrast estimate   SE  df t.ratio p.value
#F1 - F0    -27.79 3.36 567  -8.263  <.0001
#F2 - F1    -25.02 2.22 522 -11.248  <.0001
#F3 - F2      1.34 2.41 563   0.554  0.9565
#F4 - F3      7.66 3.06 548   2.506  0.0459


# Display the contrasts and their p-values
summary(pairwise_contrasts)
#Degrees-of-freedom method: kenward-roger 
#P value adjustment: mvt method for 4 tests 

# Key takeaways:
# Brood Size Declines in F2, F3, and F4: There is a significant decline in brood size in generations F2 (-25.02), F3 (-23.68), and F4 (-16.02) compared to the baseline generation (G0)
# The random effect for founder ID indicates that there is between-lineage variability in brood sizes, although the random effect variance (43.66) is much smaller than the residual variance, implying that most of the variation is within lineages rather than between them.

# Simple linear regression without random effects
#clutch_lm <- lm(brood_size ~ generation, data = sex_ratios_all)
#summary(clutch_lm)

# does inbreeding affect male or female mortality differently (i.e. is above due to lower male mortality?)
# To test if the decline in brood size differs between males and females, you can fit a model with an interaction term between generation and sex.

# GLM with binomial family to assess if generation influences male/female ratio
sex_ratios_all_long <- sex_ratios_all %>%
  pivot_longer(cols = c(males, females), 
               names_to = "sex", 
               values_to = "offspring_count")
interaction_model <- lmer(brood_size ~ generation * sex + (1 | founder_id), data = sex_ratios_all_long)

# Summary of the model
summary(interaction_model)
# ChatGPT summary: There is no significant difference in brood size between males and females, and the effect of generation on brood size does not differ by sex (as indicated by the non-significant interaction terms). This implies that the decline in brood size across generations is not sex-biased; it affects both males and females equally.

# Could plot male and female brood sizes across generations? Bar chart with error bars probably makes the most sense
# Convert 'generation' to a factor and specify the order of levels
sex_ratios_all_long$generation <- factor(sex_ratios_all_long$generation, levels = c("G0", "F1", "F2", "F3", "F4", "F5"))
# Create the boxplot
p <- ggplot(sex_ratios_all_long, aes(x = generation, y = offspring_count, fill = sex)) +
  geom_boxplot(notch = T, outlier.alpha = 0.7) +  # Create boxplot
  theme_classic() +
  labs(x = "Generation", y = "Offspring Count", fill="Sex") +
  scale_fill_manual(values = c("males" = "royalblue3", "females" = "orange")) +  # Color for male and female boxes
  theme(axis.line=element_line(linewidth=1),
        axis.ticks=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black", face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=14, face="bold"))

print(p)

ggsave(file="/Users/robertbaird/documents/writing/Ling_SR/Figure_bits/male_female_counts_over_generations.svg", plot=p, width=5.5, height=4)#

# plot total brood size
ggplot(sex_ratios_all, aes(x = generation, y = brood_size,)) +
  geom_boxplot(notch = T, outlier.alpha = 0.7) +  # Create boxplot
  labs(x = "Generation", y = "Offspring Count", 
       title = "Offspring Count by Generation") +
  theme_minimal()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

### 5. Effect of temperature on sex ratio

# To test if there is a statistically significant effect of temperature on the sex ratios of broods while weighing by brood size, a weighted linear regression model (similar to before) can be used

head(temp)
# add brood size column
temp$brood_size <- temp$males+temp$females
# add sex ratio column
temp$sex_ratio <- temp$males/temp$brood_size
temp <- temp[complete.cases(temp),]

# medians?
median(temp[which(temp$treatment=="12"),]$sex_ratio) # 0.4736842
median(temp[which(temp$treatment=="18"),]$sex_ratio) # 0.5555556
median(temp[which(temp$treatment=="25"),]$sex_ratio) # 0.5

# Fit a weighted linear regression model
model <- lm(sex_ratio ~ treatment, data = temp, weights = brood_size)

summary(model)

# temperature has almost no effect on the sex ratio.

# plot - size of points corresponds to brood size
p <- ggplot(temp, aes(x = factor(treatment), y = sex_ratio)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", width=0.5) +  # Boxplot without outliers
  geom_jitter(aes(size = brood_size), color = "black", alpha = 0.2, width = 0.1) +  # Jittered points
  theme_classic()+
  scale_size_continuous(range = c(1, 5), name = "Brood Size") +  # Adjust size scale
  labs(x = "Temperature Treatment (°C)", y = "Proportion of male offspring") +
  theme(axis.line=element_line(linewidth=1),
        axis.ticks=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black", face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=14, face="bold"))

print(p)

ggsave(file="/Users/robertbaird/documents/writing/Ling_SR/Figure_bits/temp_effect.svg", plot=p, width=5.5, height=4)#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# also:

# 1.	Check for overall mortality differences by temperature treatment:
#  Use a linear model where brood size is the dependent variable, and treatment (temperature) is the independent variable.
# 2.	Check for sex-specific mortality differences:
#  To directly test whether sex interacts with temperature in determining mortality (i.e., whether the effect of temperature on brood size differs between males and females), use a linear model with an interaction term

# 1. 
p <- ggplot(temp, aes(x = factor(treatment), y = brood_size)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Boxplot without outliers
  geom_jitter(color = "black", alpha = 0.2, width = 0.1) +  # Jittered points
  scale_size_continuous(range = c(1, 5), name = "Brood Size") +  # Adjust size scale
  labs(x = "Temperature Treatment (°C)", y = "Brood size") +
  theme_classic()+
  theme(axis.line=element_line(linewidth=1),
        axis.ticks=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black", face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=14, face="bold"))

print(p)

ggsave(file="/Users/robertbaird/documents/writing/Ling_SR/Figure_bits/temp_vs_brood_size.svg", plot=p, width=3.5, height=4)#

# Ensure treatment is a factor and set 18°C as the reference
temp$treatment <- factor(temp$treatment, levels = c("18", "12", "25"))

# Run the linear model
model <- lm(brood_size ~ treatment, data = temp)

#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  14.9143     0.7067  21.103  < 2e-16 ***
#  treatment12  -4.0933     1.0893  -3.758 0.000185 ***
#  treatment25  -0.9195     1.1466  -0.802 0.422859  

# Summarize the model
summary(model)

# mortality is higher at 12 and 25 than at 18, and it's highest at 12. the 18-25 difference is nonsignificant but the 18-12 difference is significant.

# plot effect of temperature on counts

temp_long <- temp %>%
  pivot_longer(cols = c(males, females), 
               names_to = "sex", 
               values_to = "offspring_count")
temp_long$treatment <- factor(temp_long$treatment, levels = c("12", "18", "25"))
p <- ggplot(temp_long, aes(x = sex, y = offspring_count, fill = factor(treatment))) +
  geom_boxplot(notch = T, outlier.alpha = 0.7) +  # Create boxplot
  labs(x = "Temperature treatment", y = "Offspring Count", fill="Treatment") +
  scale_fill_manual(values = c("12" = "lightblue", "18" = "lightgreen", "25" = "orange")) +  # Color for male and female boxes
  theme_classic()+
  theme(axis.line=element_line(linewidth=1),
        axis.ticks=element_line(linewidth=1),
        axis.text=element_text(size=14, colour="black", face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=14, face="bold"))

print(p)

ggsave(file="/Users/robertbaird/documents/writing/Ling_SR/Figure_bits/temp_vs_counts.svg", plot=p, width=5, height=4)#

# Linear model with interaction between temperature and sex
lm_interaction <- lm(offspring_count ~ treatment * sex, data = temp_long)

# Summarize the model
summary(lm_interaction)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            6.4921     0.4062  15.983  < 2e-16 ***
#  treatment12           -0.9855     0.6260  -1.574 0.115650    
#treatment25            0.5079     0.6590   0.771 0.440948    
#sexmales               1.9302     0.5744   3.360 0.000799 ***
#  treatment12:sexmales  -2.1223     0.8853  -2.397 0.016647 *  
#  treatment25:sexmales  -1.9353     0.9319  -2.077 0.038002 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 7.209 on 1468 degrees of freedom
#Multiple R-squared:  0.0224,	Adjusted R-squared:  0.01907 
#F-statistic: 6.727 on 5 and 1468 DF,  p-value: 3.299e-06

# Summary:
# there is no significant effect of temp on production of female offspring
# there IS a signficant effect on production of male offspring, with fewer males produced at 12, suggesting perhaps the number of males produced is more sensitive to changes in temperature







