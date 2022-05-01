import React, { useState, useEffect } from 'react';
import { useHistory } from 'react-router-dom';
import SettingsIcon from '@mui/icons-material/Settings';
import Button from '@mui/material/Button';
import ListCard from 'components/general/ListCard';
import InstitutionLeaveButton from 'components/institutions/InstitutionLeaveButton';
import styles from './institutionCard.module.css';

import getProfile from 'shared/services/profile';

function InstitutionCard({ institution, onLeft }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getProfile()
      .then(setUser);
  }, [setUser]);

  const history = useHistory();
  const navigateToInstitution = () => {
    history.push(`/sequencer/institutions/${institution.id}`);
  };

  const {
    name, avatarUrl, adminIds,
  } = institution;
  return (
    <ListCard
      title={name}
      imageURL={avatarUrl}
      enforceImage
      nextToTitle={(
        <span className={styles.accessRightIndicator}>
          {adminIds.indexOf(user.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailing={[
        adminIds.indexOf(user.id) !== -1 ? (
          <Button
            key="settings"
            endIcon={<SettingsIcon />}
            variant="outlined"
            sx={{ marginRight: '6px' }}
          >
            Settings
          </Button>
        ) : <div key="nothing" />,
        <InstitutionLeaveButton key="leave" institution={institution} onLeft={onLeft} />,
      ]}
      onClick={navigateToInstitution}
    />
  );
}

export default InstitutionCard;
