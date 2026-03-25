def validate_input(data):
    """
    Validates the input data for certain criteria.
    """
    if not isinstance(data, (int, float)):  # Check if the input is a number
        raise ValueError("Input must be a number.")
    if data < 0:  # Check for non-negative numbers
        raise ValueError("Input must be a non-negative number.")
    return True

def validate_string(input_string):
    """
    Validates that the input string is not empty and is alphanumeric.
    """
    if not input_string or not isinstance(input_string, str):
        raise ValueError("Input must be a non-empty string.")
    if not input_string.isalnum():
        raise ValueError("Input must be alphanumeric.")
    return True

# Additional validation functions can be added here.