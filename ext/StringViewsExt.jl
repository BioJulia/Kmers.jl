module StringViewsExt

using StringViews: StringView
using Kmers: Kmers

# This extension is important because FASTX uses string views.
# The documentation of StringViews promises that the underlying
# string is UTF-8 encoded.
Kmers.is_ascii(::Type{<:StringView}) = true

end # module
