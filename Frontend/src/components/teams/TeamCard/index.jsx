import React from 'react';
import { useHistory } from 'react-router-dom';
import SettingsIcon from '@mui/icons-material/Settings';
import ListCard from 'components/general/ListCard';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';
import Button from 'components/CustomButton';
import styles from './teamCard.module.css';
import { useAuth } from 'shared/context/authContext';

function TeamCard({ team, onLeft }) {
  const [user] = useAuth();

  const history = useHistory();
  const navigateToTeam = () => history.push(`/sequencer/teams/${team.id}`);

  const preventBubble = (evt) => {
    evt.stopPropagation();
  };

  const {
    name, description, adminIds, memberIds,
  } = team;
  return (
    <ListCard
      title={name}
      description={description}
      nextToTitle={(
        <span className={styles.accessRightIndicator}>
          {adminIds.indexOf(user._id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailing={[
        adminIds.indexOf(user._id) !== -1 ? (
          <Button
            key="settings"
            endIcon={<SettingsIcon />}
            variant="outlined"
            sx={{ marginRight: '6px' }}
            onClick={navigateToTeam}
          >
            Settings
          </Button>
        ) : <div key="nothing" />,
        memberIds?.includes(user._id) && (
        <div key="leave" onClick={preventBubble} onKeyPress={preventBubble} role="button" tabIndex={0}>
          <TeamLeaveButton team={team} onLeft={onLeft} />
        </div>
        ),
      ]}
      onClick={navigateToTeam}
    />
  );
}

export default TeamCard;
