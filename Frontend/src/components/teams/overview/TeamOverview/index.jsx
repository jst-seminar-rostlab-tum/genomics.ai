import React, { useState, useEffect } from 'react';
import { Alert } from '@mui/material';
import Stack from '@mui/material/Stack';
import Typography from '@mui/material/Typography';
import PlusIcon from 'components/general/PlusIcon';
import TeamCard from 'components/teams/TeamCard';
import styles from './teamOverview.module.css';
import TeamService from 'shared/services/Team.service';
import TeamCreationDialog from 'components/teams/overview/TeamCreationDialog';

function TeamOverview() {
  const [teams, setTeams] = useState([]);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    setIsLoading(true);
    TeamService.getMyTeams()
      .then(setTeams)
      .finally(() => setIsLoading(false));
  }, []);

  function onLeft(team) {
    setTeams(teams.filter((i) => i.id !== team.id));
  }

  const [createOpen, setCreateOpen] = useState(false);

  return (
    <>
      <Stack direction="row" className="stack" alignItems="Center">
        <Typography variant="h5" sx={{ pr: 1 }}>Your Teams</Typography>
        <PlusIcon onClick={() => setCreateOpen(true)} />
      </Stack>
      <br />
      {!isLoading && teams.length === 0 ? (
        <Alert severity="info">
          You have not created any teams yet. Create one by clicking the Plus-Icon.
        </Alert>
      ) : ''}
      {teams.map((team) => (
        <div key={team.id}>
          <TeamCard
            team={team}
            onLeft={(t) => onLeft(t)}
          />
          <div className={styles.cardSpacing} />
        </div>
      ))}
      {createOpen && <TeamCreationDialog
        open={createOpen}
        handleClose={() => setCreateOpen(false)}
        onCreated={(newTeam) => setTeams([...teams, newTeam])}
      /> }
    </>
  );
}

export default TeamOverview;
