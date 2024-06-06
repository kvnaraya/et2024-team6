
#read in data
rm(list=ls(all=TRUE))

##Update packages
update.packages()

##Load relevant package_cs
library(reshape2)
library(nlme)    #Non-linear mixed effects package
library(car)   #Regression package
library(lme4)    #Linear mixed effects package
library(ggplot2)   #Advanced plotting package
library(language_cR)  ##Package_cs containing useful functions for language research 
library(lattice)  ##Plotting package
library(lmerTest) 
library(dplyr)
library(tidyr)
library(esquisse)
library(readxl)
library(ggplot2)
library(extrafont)
setwd("/Users/keertananarayanan/Desktop/et2024")
data = read.table("timecourse.csv", header = TRUE, sep = ",")
data1 <- read.table("timecourse.csv", header = TRUE, sep =",")
data1$wordtype <- factor(data1$wordtype, levels = c("target", "competitor", "unrelated"))
data2 = read_excel("timecourse2.xlsx")
data4 = read_excel("timecourse3.xlsx")
data3 <- read.table("timecoursebyage.csv", header = TRUE, sep = ",")
data3$agegroup <- factor(data3$agegroup, levels = c("(7-8)", "(11-12)", "(16-17)"))

#Register the fonts for use
# Import the Adobe Type 1 font
font_import(pattern = "Times", prompt = FALSE)
# Load the font for use
loadfonts()

# reorder 
data$wordtype <- factor(data$wordtype, levels = c("target", "competitor", "unrelated"))

library(ggplot2)
library(grid)

plot1 <- ggplot(data) +
  aes(x = Time, y = value, group = wordtype) +
  geom_line(aes(linetype = wordtype), linewidth = 1.5) + 
  geom_vline(xintercept = 300, linetype = "solid", color = "black", linewidth = 0.5) +  # Add the vertical line
  scale_linetype_discrete(name = "Word Type", labels = function(x) stringr::str_to_title(x)) +  # Capitalize first letter
  labs(x = "Time (msec)", y = "Proportion of Fixations") +
  theme_minimal() +
  facet_wrap(vars(Trialcodetype)) +
  theme(legend.justification = c(1, 0.5), legend.position = c(1, 0.5),
        strip.text.x = element_text(size = 18, family = "Times New Roman", color = "black"),
        axis.title.x = element_text(size = 18, family = "Times New Roman", color = "black"),
        axis.text.x = element_text(size = 18, family = "Times New Roman", color = "black"),
        axis.title.y = element_text(size = 18, family = "Times New Roman", color = "black"),
        axis.text.y = element_text(size = 18, family = "Times New Roman", color = "black"),
        legend.title = element_text(size = 18, family = "Times New Roman", color = "black"),
        legend.text = element_text(size = 18, family = "Times New Roman", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5)) + 
  scale_y_continuous(limits = c(0, 1.0), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 2000), expand = c(0, 0)) + 
  theme(plot.margin = margin(l = 0.2, r = 0.6, b = 0.5, t = 0.5, unit = "cm"),
        plot.background = element_rect(fill = "white"),  # Adding white background to plot
        panel.spacing = unit(2.0, "lines")) +  # Adjusting panel spacing
  guides(linetype = guide_legend(override.aes = list(size = 1.5)))  # Adjust legend linetype size

print(plot1)

# Save the plot as an EPS file with Times New Roman font
ggsave("Figure1.pdf", plot = plot1, width = 10, height = 6, units = "in")

# Figure 2
#SLACHL by age group and just cohort 

cohort <- subset(data3, Trialcodetype== "Cohort")

figure2 <- ggplot(cohort) +
  aes(x = Time, y = fixations, colour = agegroup) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(name = "Age Group", values = c(`(7-8)` = "#1204E9", `(11-12)` = "#E90606", `(16-17)` = "#F48408")) +
  labs(x = "Time (msec)", y = "Proportion of Fixations") + 
  theme_minimal() +
  facet_wrap(vars(type), scales = "free_y") + 
  theme(legend.justification = c(1, 0.5), legend.position = c(1, 0.5),
        strip.text.x = element_text(size = 18, family = "Times New Roman", color = "black"),
        axis.title.x = element_text(size = 18, family = "Times New Roman", color = "black"),
        axis.text.x = element_text(size = 18, family = "Times New Roman", color = "black"),
        axis.title.y = element_text(size = 18, family = "Times New Roman", color = "black"),
        axis.text.y = element_text(size = 18, family = "Times New Roman", color = "black"),
        legend.title = element_text(size = 18, family = "Times New Roman", color = "black"),
        legend.text = element_text(size = 18, family = "Times New Roman", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.line.x = element_blank(), 
        axis.ticks = element_line(color = "black", linewidth = 0.5)) +
  geom_hline(yintercept = 0) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 2000), expand = c(0, 0)) +
  theme(plot.margin = margin(l= .2, r = .6, unit = "cm"))  # Add 1 centimeter of white space below and to the right of the plot
# Save the plot as an EPS file with Times New Roman font
ggsave("Figure2.pdf", figure2, device = "pdf", family = "Times New Roman", width = 10, height = 6)

# Figure 3 
# VVWP total
library(ggplot2)
library(patchwork)
data4$type <- factor(data4$type, levels = c("Target", "Competitor", "Unrelated"))
# Create the first plot
plot1 <- ggplot(data4, aes(x = Time, y = fixations, group = type)) +
  geom_line(aes(linetype = type), linewidth = 1.5) + 
  scale_linetype_manual(name = "Word Type", values = c(Target = "solid", Competitor = "dashed", Unrelated = "22" )) +
  labs(y = "Proportion of Fixations") +
  theme_minimal() +
  facet_wrap(vars(type2), scales = "free") + 
  theme(legend.justification = c(1, .5), legend.position = c(1, .5),
        strip.text.x = element_text(size = 24, family = "Times New Roman", color = "black"),  # Adjusted font size and color
        axis.title.x = element_blank(),  # Adjusted font size
        axis.text.x = element_text(size = 24, family = "Times New Roman", margin = margin(t = 0.2, unit = "cm"), color = "black"),  # Adjusted font size and color
        axis.title.y = element_text(size = 24, family = "Times New Roman", color = "black"),  # Adjusted font size and color
        axis.text.y = element_text(size = 24, family = "Times New Roman", color = "black"),   # Adjusted font size and color
        legend.title = element_text(size = 24, family = "Times New Roman", color = "black"),  # Adjusted font size and color
        legend.text = element_text(size = 24, family = "Times New Roman", color = "black"),   # Adjusted font size and color
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(),
        axis.ticks = element_line(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(linewidth = 0.5)) +  # Show x-axis ticks
  geom_hline(yintercept = 0) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 2000, by = 500)) +
  expand_limits(x = c(0, 2000)) +  # Ensure no x-axis text is cut off
  theme(plot.margin = margin(b = 1, unit = "cm"))  # Add 1 centimeter of white space below the plot

# Create the second plot
plot2 <- ggplot(data2, aes(x = Time, y = fixations, colour = type)) +
  geom_line(linewidth = 1.5) + 
  scale_color_manual(name = "Age Group", values = c('(7-8)' = "#0F3BE4", '(11-12)' = "#E60202", '(16-17)' = "#ED7306")) + 
  labs(y = "Proportion of Fixations") +
  theme_minimal() +
  facet_wrap(vars(type2), scales = "free") + 
  theme(legend.justification = c(1, .5), legend.position = c(1, .5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 24, family = "Times New Roman", color = "black"),   # Adjusted font size and color
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 24, family = "Times New Roman", color = "black"),   # Adjusted font size and color
        legend.title = element_text(size = 24, family = "Times New Roman", color = "black"),  # Adjusted font size and color
        legend.text = element_text(size = 24, family = "Times New Roman", color = "black"),   # Adjusted font size and color
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(),
        axis.ticks = element_line(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_line(linewidth = 0.5),  # Show x-axis ticks
        strip.text = element_text(size = 24, family = "Times New Roman", color = "black")) +  # Adjusted font size and color
  geom_hline(yintercept = 0) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 2000, by = 500)) + 
  expand_limits(x = c(0, 2000)) +  # Ensure no x-axis text is cut off
  theme(plot.margin = margin(b = 1, r = 1, unit = "cm"))  # Add 1 centimeter of white space below and to the right of the plot

# Combine the plots side by side
combined_plot <- plot1 + plot2 + plot_layout(widths = c(.5, 1), nrow = 1)

# Save the combined plot with annotation
ggsave("CombinedPlots.pdf", combined_plot, device = "pdf", width = 20, height = 8, pointsize = 30)

# Figure 6
semantic <- subset(data3, Trialcodetype== "Semantic")
semantic <- subset(semantic, type== "B. Competitor Effect")

figure6 <- ggplot(semantic) +
  aes(x = Time, y = fixations, colour = agegroup) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(name = "Age Group", values = c(`(7-8)` = "#1204E9", `(11-12)` = "#E90606", `(16-17)` = "#F48408")) +
  labs(x = "Time (msec)", y = "Proportion of Fixations") + 
  theme_minimal() +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1),  # Move legend to top right
        strip.text.x = element_text(size = 18, family = "Times New Roman", color = "black"),  # Set strip text color to black
        axis.title.x = element_text(size = 18, family = "Times New Roman", color = "black"),  # Set axis title color to black
        axis.text.x = element_text(size = 18, family = "Times New Roman", color = "black"),   # Set x-axis text color to black
        axis.title.y = element_text(size = 18, family = "Times New Roman", color = "black"),  # Set axis title color to black
        axis.text.y = element_text(size = 18, family = "Times New Roman", color = "black"),   # Set y-axis text color to black
        legend.title = element_text(size = 18, family = "Times New Roman", color = "black"),  # Set legend title color to black
        legend.text = element_text(size = 18, family = "Times New Roman", color = "black"),   # Set legend text color to black
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"),  # Set axis line color to black
        axis.ticks = element_line(color = "black", linewidth = 0.5),  # Set axis ticks color and size
        axis.ticks.length = unit(0.2, "cm"),  # Set tick length
        axis.line.x = element_line(),  # Show x-axis line
        strip.text = element_text(size = 18, family = "Times New Roman", color = "black")) +  # Set strip text color
  geom_hline(yintercept = 0) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, 2000), expand = c(0, 0), breaks = seq(0, 2000, by = 500), minor_breaks = seq(0, 2000, by = 100)) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5)) + # Center legend title 
  theme(plot.margin = margin(l = .4, r = .8, unit = "cm")) 

ggsave("Figure6.pdf", figure6, device = "pdf", family = "Times New Roman", width = 10, height = 6)




