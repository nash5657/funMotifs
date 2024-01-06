# Step 1: Install Java
# Skip this step if Java 8 or above is already installed
# Example for installing OpenJDK 8
sudo apt update
sudo apt install openjdk-8-jdk

# Step 2: Install Scala
# Skip this step if Scala is already installed
# Example for installing Scala
sudo apt install scala

# Step 3: Download and Install Apache Spark
# Replace the version numbers with the latest ones if needed
SPARK_VERSION="3.2.0"
HADOOP_VERSION="3.2"
# Download Spark
wget https://downloads.apache.org/spark/spark-${SPARK_VERSION}/spark-${SPARK_VERSION}-bin-hadoop${HADOOP_VERSION}.tgz
# Unpack the Spark archive
tar -xzf spark-${SPARK_VERSION}-bin-hadoop${HADOOP_VERSION}.tgz
# Move Spark to /opt
sudo mv spark-${SPARK_VERSION}-bin-hadoop${HADOOP_VERSION} /opt/spark
# Create a symbolic link
sudo ln -s /opt/spark /opt/spark-3.2.0

# Step 4: Set up environment variables
# Add the following lines to the end of your ~/.bashrc or ~/.zshrc file
echo 'export SPARK_HOME=/opt/spark' >> ~/.bashrc
echo 'export PATH=$PATH:$SPARK_HOME/bin' >> ~/.bashrc
# Load the changes
source ~/.bashrc

# Step 5: Verify the installation
# Start the Spark shell
spark-shell
