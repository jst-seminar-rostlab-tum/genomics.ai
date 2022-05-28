import React, { useEffect, useState } from 'react';
import InstitutionService from 'shared/services/Institution.service';

const InstitutionContext = React.createContext();
InstitutionContext.displayName = 'Institution Context';

function InstitutionProvider(props) {
  const [institutions, setInstitutions] = useState([]);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    setIsLoading(true);
    InstitutionService.getMyInstitutions()
      .then((newInstitutions) => {
        setInstitutions(newInstitutions);
        setIsLoading(false);
      })
      .finally(() => setIsLoading(false));
  }, []);

  const value = { isLoading, institutions, setInstitutions };

  return (
    // eslint-disable-next-line react/jsx-props-no-spreading
    <InstitutionContext.Provider value={value} {...props} />
  );
}

function useInstitutions() {
  const context = React.useContext(InstitutionContext);
  if (context === undefined) {
    throw new Error('useInstitutions must be used within a InstitutionContext');
  }
  return context;
}

export { InstitutionProvider, useInstitutions };
