import logging
from typing import Any, Dict, Union

# Set up logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def validate_input(data: Dict[str, Any]) -> bool:
    """
    Validates the input data.
    
    :param data: Input data to validate.
    :return: True if valid, False otherwise.
    """
    if not isinstance(data, dict):
        logger.error("Input data is not a dictionary.")
        return False
    # Add more validation logic as necessary.
    logger.info("Input data is valid.")
    return True


def process_data(data: Dict[str, Any]) -> Union[None, str]:
    """
    Processes the input data and returns a string result.
    
    :param data: Input data to process.
    :return: Result as a string if successful, None otherwise.
    """
    if not validate_input(data):
        logger.error("Invalid input data.")
        return None
    # Processing logic here.
    result = 'Processed data'  # Placeholder for real processing logic.
    logger.info("Successfully processed the data.")
    return result
