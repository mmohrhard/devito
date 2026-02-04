#ifndef _DEVITO_CUDA_LOGGING_H
#define _DEVITO_CUDA_LOGGING_H

#include <string>
#include <memory>
#include <stdexcept>


#ifndef STRINGIFY
#define STRINGIFY(x) _stringify(x)
#define _stringify(x) #x
#endif

/**
 * LOGGING
 */

enum LogLevel {
  CRITICAL = 50,
  FATAL = 50,
  ERROR = 40,
  WARNING = 30,
  WARN = 30,
  INFO = 20,
  DEBUG = 10,
  NOTSET = 0
};

typedef void (*LogHandler)(int logLevel, const char *message);

LogHandler _logHandler = nullptr;
LogLevel _minLogLevel = LogLevel::DEBUG;

extern "C" void setLogHandler(void *handler) {
  _logHandler = (LogHandler)handler;
}

extern "C" void setLogLevel(int logLevel) {
    _minLogLevel = (LogLevel)logLevel;
}

template <typename... Args>
std::string string_format(const std::string &format, Args... args) {
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) +
               1; // Extra space for '\0'
  if (size_s <= 0) {
    throw std::runtime_error("Error during formatting.");
  }
  auto size = static_cast<size_t>(size_s);
  std::unique_ptr<char[]> buf(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args...);
  return std::string(buf.get(),
                     buf.get() + size - 1); // We don't want the '\0' inside
}

template <typename... Args>
inline void log(int logLevel, const std::string &format, Args... args) {
  if (_logHandler == nullptr || logLevel < _minLogLevel)
    return;

  std::string message = string_format(format, std::forward<Args>(args)...);

  _logHandler(logLevel, message.c_str());
}

template <typename... Args>
inline void debug(const std::string &format, Args... args) {
  log(LogLevel::DEBUG, format, std::forward<Args>(args)...);
}

template <typename... Args>
inline void info(const std::string &format, Args... args) {
  log(LogLevel::INFO, format, std::forward<Args>(args)...);
}

template <typename... Args>
inline void warn(const std::string &format, Args... args) {
  log(LogLevel::WARNING, format, std::forward<Args>(args)...);
}

template <typename... Args>
inline void critical(const std::string &format, Args... args) {
  log(LogLevel::CRITICAL, format, std::forward<Args>(args)...);
}

#ifdef DEVITO_CUDA_VERBOSE_MPI_REGIONS
#define LOG_COMPUTE(NAME, ...) { \
    debug("compute %s", STRINGIFY(NAME)); \
}

#define LOG_REMAINDER(NAME, NREGIONS, ...) { \
    debug("remainder %s with %d regions", STRINGIFY(NAME), NREGIONS); \
}
#else
#define LOG_COMPUTE(...) {}
#define LOG_REMAINDER(...) {}
#endif

#endif // _DEVITO_CUDA_LOGGING_H
