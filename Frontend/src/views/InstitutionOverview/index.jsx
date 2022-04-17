import React, { useState, useEffect } from 'react';
import HeaderView from 'components/HeaderView';
import InstitutionCard from 'components/InstitutionCard';
import styles from './institutionOverview.module.css';
import queryMyInstitutions from 'shared/services/mock/institutions';

function InstitutionOverview({ sidebarShown }) {
  const [institutions, setInstitutions] = useState([]);

  useEffect(() => {
    const updateInstitutions = () => queryMyInstitutions()
      .then((newInstitutions) => setInstitutions(newInstitutions))
      .catch((ignored) => { console.log(ignored); });
    updateInstitutions();
  }, [setInstitutions]);
  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title="My Institutions"
    >
      {institutions.map((institution) => (
        <InstitutionCard key={institution.id} institution={institution} />
      ))}
    </HeaderView>
  );
}

export default InstitutionOverview;
