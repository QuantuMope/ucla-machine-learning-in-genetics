Use shell script **phaser** to run.

Takes following commands:
- impute $(masked_data) - imputes masked data and outputs to text file called output
- compare $(true_data) $(output) - outputs different num of different alleles between output and true data
- total_mask $(masked_data) - outputs the num of masked alleles
- run $(output) - runs the phaser on imputed data