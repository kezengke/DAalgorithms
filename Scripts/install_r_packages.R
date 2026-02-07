# Install required R packages for differential abundance analysis

# List of required packages
required_packages <- c(
  "DESeq2",
  "edgeR", 
  "ALDEx2",
  "ANCOMBC",
  "metagenomeSeq",
  "BiocManager"
)

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "edgeR", "ALDEx2", "ANCOMBC", "metagenomeSeq"), 
                     update = FALSE, ask = FALSE)

# Install CRAN packages if needed
if (!require("ALDEx2", quietly = TRUE)) {
  install.packages("ALDEx2")
}

# Check which packages are installed
cat("Checking installed packages:\n")
for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("✓", pkg, "is installed\n"))
  } else {
    cat(paste("✗", pkg, "is NOT installed\n"))
  }
}

