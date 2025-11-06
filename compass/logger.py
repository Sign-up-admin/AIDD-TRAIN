import os
import datetime


class TrainingLogger:
    """
    A simple logger for training processes.

    Handles logging of general messages, warnings, and errors to both the console
    and dedicated log files.
    """

    def __init__(self, log_dir="logs"):
        """
        Initializes the logger.

        Creates the log directory and sets up the log files.

        Args:
            log_dir (str): The directory where log files will be stored.
        """
        os.makedirs(log_dir, exist_ok=True)
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

        # Main log file for all messages
        self.log_file = os.path.join(log_dir, f"training_log_{timestamp}.txt")

        # Dedicated log file for warnings and errors
        self.error_log_file = os.path.join(log_dir, f"training_errors_{timestamp}.txt")

        self.log(f"Main log file created at {self.log_file}")
        self.log(f"Error log will be saved to {self.error_log_file}")

    def log(self, message):
        """
        Appends a general, timestamped message to the main log file.

        Args:
            message (str): The message to log.
        """
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_message = f"[{timestamp}] [INFO] {message}\n"

        print(message)

        with open(self.log_file, "a", encoding="utf-8") as f:
            f.write(log_message)

    def log_warning(self, message):
        """
        Appends a timestamped warning message to both the main and error log files.

        Args:
            message (str): The warning message to log.
        """
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        warning_message = f"[{timestamp}] [WARNING] {message}\n"

        print(f"WARNING: {message}")

        with open(self.log_file, "a", encoding="utf-8") as f:
            f.write(warning_message)
        with open(self.error_log_file, "a", encoding="utf-8") as f:
            f.write(warning_message)

    def log_error(self, message):
        """
        Appends a timestamped error message to both the main and error log files.

        Args:
            message (str): The error message to log.
        """
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        error_message = f"[{timestamp}] [ERROR] {message}\n"

        print(f"ERROR: {message}")

        with open(self.log_file, "a", encoding="utf-8") as f:
            f.write(error_message)
        with open(self.error_log_file, "a", encoding="utf-8") as f:
            f.write(error_message)
