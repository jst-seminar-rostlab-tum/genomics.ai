import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import HeaderView from 'components/general/HeaderView';
import TeamProjectList from 'components/teams/detail/TeamProjectList';
import TeamMemberList from 'components/teams/detail/TeamMemberList';
import TeamAdminHeaderRight from 'components/teams/detail/TeamAdminHeaderRight';
import TeamUserHeaderRight from 'components/teams/detail/TeamUserHeaderRight';
import TeamHeaderOptions from 'components/teams/detail/TeamHeaderOptions';
import TeamService from 'shared/services/Team.service';
import InstitutionService from 'shared/services/Institution.service';
import TeamInviteButton from 'components/teams/detail/TeamInviteButton';
import TextField from '@mui/material/TextField';
import { useAuth } from 'shared/context/authContext';
import Button from 'components/CustomButton';

export default function TeamPage() {
  const { id } = useParams();
  const [team, setTeam] = useState({});
  const [user] = useAuth();
  const [institution, setInstitution] = useState({});
  const [newDescription, setNewDescription] = useState("");
  const [descriptionChanged, setDescriptionChanged] = useState(false);

  const handleDescriptionChange = (event) => {
    event.preventDefault();
    setDescriptionChanged(true);
    setNewDescription(event.target.value);
  };

  function updateTeam() {
    TeamService.getTeam(id)
      .then(setTeam)
      .catch((ignored) => { console.error(ignored); });
  };

  async function updateVisibility(newVisibility) {
    await TeamService.changeTeamVisibility(id, newVisibility);
    updateTeam()
  }

  async function updateDescription(newDescription) {
    await TeamService.changeTeamDescription(id, newDescription);
    updateTeam()
    setDescriptionChanged(false);
  }

  useEffect(() => {
    if (id == null) return;
    updateTeam()
  }, [setTeam]);

  useEffect(() => {
    setNewDescription(team.description);
  }, [team])

  // Institution may be undefined
  useEffect(() => {
    InstitutionService.getInstitution(team.institutionId)
      .then((newInstitution) => setInstitution(newInstitution))
      .catch((ignored) => { console.error(ignored); });
  }, [team, setInstitution]);

  const isAdmin = team.adminIds ? team.adminIds.includes(user._id) : false;

  if (!team) return <></>;

  return (
    <HeaderView
      title={team.name}
      rightOfTitle={(
        <TeamHeaderOptions
          team={team}
          isAdmin={isAdmin}
          institution={institution}
        />
      )}
      replaceHeaderRight={
        (isAdmin && <TeamAdminHeaderRight team={team} updateVisibility={updateVisibility} />)
        || <TeamUserHeaderRight institution={institution} team={team} user={user} updateTeam={updateTeam} />
      }
    >
      <br />
      <section>
        <h2>Description</h2>
        <hr />
        <div>
          <TextField
            id="description"
            multiline
            minRows={3}
            maxRows={5}
            value={newDescription}
            InputProps={{
              readOnly: !isAdmin,
            }}
            fullWidth
            onChange={handleDescriptionChange}
            variant="standard"
          />
          {descriptionChanged
            && (
              <Button variant="outlined" type="secondary" onClick={() => updateDescription(newDescription)}>
                Submit
              </Button>
            )}
        </div>
      </section>
      <section>
        <h2>Projects</h2>
        <hr />
        <TeamProjectList teamId={id} />
      </section>
      <section>
        <h2>Members</h2>
        <hr />
        <TeamMemberList
          team={team}
          onMemberRemoved={(_team, member) => {
            setTeam({
              ...team,
              memberIds: team.memberIds.filter((mId) => mId !== member.id),
              adminIds: team.adminIds.filter((aId) => aId !== member.id),
            });
          }}
          onMakeAdmin={(_team, member) => {
            setTeam({
              ...team,
              memberIds: team.memberIds.filter((mId) => mId !== member.id),
              adminIds: [...team.adminIds, member.id],
            });
          }}
          onRemoveAdmin={(_team, member) => {
            setTeam({
              ...team,
              memberIds: [...team.memberIds, member.id],
              adminIds: team.adminIds.filter((aId) => aId !== member.id),
            });
          }}
        />
      </section>
      {isAdmin && <TeamInviteButton team={team} />}
    </HeaderView>
  );
}
