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
      <div className={styles.content}>
        {institutions.map((institution) => (
          <div key={institution.id}>
            <InstitutionCard institution={institution} />
            <div className={styles.cardSpacing} />
          </div>
        ))}
      </div>
    </HeaderView>
  );
}

export default InstitutionOverview;
