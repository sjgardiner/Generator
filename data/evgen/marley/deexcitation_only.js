{
  // Used only by genie::MarleyDeExcitation via genie::MarleyInterface/DeExcitationOnly.
  // No primary interaction is generated here -- GENIE (interaction + INCL FSI)
  // supplies the excited remnant nucleus per-event; MARLEY only de-excites it.
  // Explicit null (not omitted) so MARLEY logs and skips this setup instead
  // of erroring, matching the "Null reactions array.../Null source
  // specification..." behavior described for this interface.
  reactions: null,
  source: null,
}

