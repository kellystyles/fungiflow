import argparse

# Function to parse and validate command line arguments

def parse_arguments():
    parser = argparse.ArgumentParser(description='Argument parsing with type validation.')

    # Add your arguments here with type validation
    parser.add_argument('--arg1', type=int, required=True, help='An integer argument.')
    parser.add_argument('--arg2', type=str, required=True, help='A string argument.')
    parser.add_argument('--arg3', type=float, required=False, default=0.0, help='A float argument (default: 0.0).')

    args = parser.parse_args()

    # Ensure the arguments meet any additional criteria
    validate_arguments(args)

    return args

# Function to validate the arguments after parsing

def validate_arguments(args):
    # Add your validation logic here
    if args.arg1 < 0:
        raise ValueError('arg1 must be a non-negative integer.')

if __name__ == '__main__':
    args = parse_arguments()