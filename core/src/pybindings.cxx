#include <pybindings.h>
#include <container_pybindings.h>

#include <string>
#include <exception>
#ifdef __FreeBSD__
#include <sys/endian.h>
#endif

std::string check_buffer_format(std::string fmt) {
	// Consume endian definition
	const char *format = &fmt[0];
	if (format[0] == '@' || format[0] == '=')
		format++;
#if BYTE_ORDER == LITTLE_ENDIAN
	else if (format[0] == '<')
		format++;
	else if (format[0] == '>' || format[0] == '!')
		throw std::runtime_error("Does not support big-endian numpy arrays");
#else
	else if (format[0] == '<')
		throw std::runtime_error("Does not support little-endian numpy arrays");
	else if (format[0] == '>' || format[0] == '!')
		format++;
#endif

	return std::string(format);
}

// The following implements the headerless module registration code
typedef std::map<std::string, std::deque<module_reg_func_t> > module_reg_t;
static std::unique_ptr<module_reg_t> modregs;

G3ModuleRegistrator::G3ModuleRegistrator(const char *mod, module_reg_func_t reg)
{
	if (!modregs)
		modregs = std::unique_ptr<module_reg_t>(new module_reg_t);
	log_debug("Adding registrar for module %s", mod);
	(*modregs)[mod].push_back(reg);
}

void G3ModuleRegistrator::CallRegistrarsFor(const char *mod, py::module_ &scope)
{
	auto &regs = (*modregs)[mod];

	// Call registered functions only once
	while (!regs.empty()) {
		auto reg = regs.front();
		regs.pop_front();
		reg(scope);
	}
}
