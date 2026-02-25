import matplotlib.pyplot as plt

# Parameters
file_path = "your_file.txt"
block_size = 200
line_index = 4         # 5th entry (0-based index)
column_index = 1       # e.g., second column (0-based)

# Read all lines
with open(file_path) as f:
    lines = f.readlines()

# Extract desired values
values = []
for i in range(line_index, len(lines), block_size):
    columns = lines[i].strip().split()  # Use .split(',') if it's comma-separated
    value = float(columns[column_index])  # Change float() to str/int if needed
    values.append(value)

# Plot the extracted values
plt.plot(values, marker='o')
plt.xlabel("Block Number")
plt.ylabel(f"Column {column_index + 1} Value")
plt.title("Selected Column from 5th Line of Each 200-Line Block")
plt.grid(True)
