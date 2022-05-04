import React from 'react';
import { useHistory } from 'react-router-dom';
import SettingsIcon from '@mui/icons-material/Settings';
import Button from '@mui/material/Button';
import ListCard from 'components/general/ListCard';
import InstitutionLeaveButton from 'components/institutions/InstitutionLeaveButton';
import styles from './institutionCard.module.css';

import { useAuth } from 'shared/context/authContext';

function InstitutionCard({
  institution, onLeft, replaceTrailing, disableLink,
}) {
  const [user] = useAuth();

  const history = useHistory();
  const navigateToInstitution = () => {
    if (disableLink) return;
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
      trailing={replaceTrailing ? replaceTrailing(institution) : [
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
