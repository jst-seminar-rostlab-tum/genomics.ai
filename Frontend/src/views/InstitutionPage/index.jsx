import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import TextField from '@mui/material/TextField';
import EditIcon from '@mui/icons-material/ModeEditOutline';
import SaveOutlinedIcon from '@mui/icons-material/SaveOutlined';
import InstitutionMemberList from 'components/institutions/InstitutionMemberList';
import styles from './institutionPage.module.css';
import CircularProgress from '@mui/material/CircularProgress';
import InstitutionTeamList from 'components/institutions/InstitutionTeamList';
import InstitutionAvatar from 'components/institutions/InstitutionAvatar';
import InstitutionBackgroundImageUploadDialog from 'components/general/upload/InstitutionBackgroundImageUploadDialog';
import { useAuth } from 'shared/context/authContext';
import InstitutionService from 'shared/services/Institution.service';
import defaultBackgroundPicture from 'assets/institution-default-background.jpg';
import InstitutionInviteButton from 'components/institutions/InstitutionInviteButton';

function InstitutionPage() {
  const { id } = useParams();
  const [institution, setInstitution] = useState({});
  const [user] = useAuth();
  const [institutionLoaded, setInstitutionLoaded] = useState(false);
  const [backgroundUploadOpen, setBackgroundUploadOpen] = useState(false);

  function isAdmin() {
    return (institution.adminIds || []).includes(user._id);
  }

  function updateInstitution() {
    InstitutionService.getInstitution(id)
      .then(setInstitution)
      .catch((ignored) => { console.error(ignored); });
  }

  function editDetails() {
    InstitutionService.updateDetails(id, institution.description);
  }

  const handleDescriptionChange = (event) => {
    setInstitution({
      ...institution,
      description: event.target.value,
    });
  };

  useEffect(() => {
    InstitutionService.getInstitution(id)
      .then((newInstitution) => {
        setInstitution(newInstitution);
        setInstitutionLoaded(true);
      });
  }, []);

  function onLeft(/* team */) {
    // setTeams(teams.filter((i) => i.id !== team.id));
  }

  if (!institutionLoaded) {
    return <CircularProgress />;
  }

  return (
    <>
      <div
        className={styles.background}
        style={{
          backgroundImage: `url(${institution.backgroundPictureURL || defaultBackgroundPicture})`,
          resizeMode: 'stretch',
        }}
      >
        <div className={styles.institutionIcon}>
          <InstitutionAvatar
            institution={institution}
            editable={isAdmin()}
            onChange={(newUrl) => {
              // update without reload
              console.log("change", newUrl, institution);
              setInstitution({ ...institution, profilePictureURL: newUrl });
            }}
          />
        </div>
        <div
          style={{
            display: 'flex',
            justifyContent: 'center',
            marginBottom: '10px',
            marginTop: '10px',
          }}
        >
          <TextField
            minRows={1}
            maxRows={1}
            value={institution.name}
            sx={{
              input: {
                textAlign: 'center',
                color: 'white',
                fontSize: 40,
                backgroundColor: 'rgba(0,38,68,0.5)',
                backdropFilter: 'blur(10px)',
                borderRadius: '34px',
              },
            }}
            InputProps={{
              disableUnderline: true,
              readOnly: true,
            }}
            style={{ width: '700px' }}
            variant="standard"
          />
        </div>
        <div style={{ display: 'flex', justifyContent: 'center' }}>
          <TextField
            minRows={1}
            maxRows={1}
            value={institution.country}
            sx={{
              input: {
                textAlign: 'center',
                color: 'white',
                fontSize: 25,
                backgroundColor: 'rgba(0,38,68,0.5)',
                backdropFilter: 'blur(10px)',
                borderRadius: '23px',
              },
            }}
            InputProps={{
              disableUnderline: true,
              readOnly: true,
            }}
            style={{ width: '300px' }}
            variant="standard"
          />
        </div>
        <p className={styles.imageText}>
          <span>
            {institution.memberIds?.length}
            {' Members'}
          </span>
        </p>
        {isAdmin() && (
          <button
            className={styles.bgImgEditButton}
            type="button"
            onClick={() => setBackgroundUploadOpen(true)}
          >
            <span>Edit Background</span>
            <EditIcon fontSize="small" />
          </button>
        )}
      </div>
      <div className={styles.test}>
        <section>
          {isAdmin() && (
            <button
              className={styles.editDetailsButton}
              type="button"
              onClick={() => editDetails()}
            >
              <span>Save Edits</span>
              <SaveOutlinedIcon fontSize="small" />
            </button>
          )}
          <h2>Description</h2>
          <hr />
          <TextField
            multiline
            minRows={3}
            maxRows={5}
            value={institution.description}
            InputProps={{
              readOnly: !isAdmin(),
            }}
            style={{ width: '100%' }}
            onChange={handleDescriptionChange}
            variant="standard"
          />
        </section>
        <section>
          <h2>Teams</h2>
          <hr />
          <div className={styles.content}>
            <InstitutionTeamList
              onLeft={(t) => onLeft(t)}
              institution={institution}
            />
            <div className={styles.cardSpacing} />
          </div>
        </section>
        <section>
          <h2>Members</h2>
          <hr />
          <InstitutionMemberList
            institution={institution}
            // eslint-disable-next-line react/jsx-no-bind
            updateInstitution={updateInstitution}
          />
        </section>
      </div>
      {isAdmin() && <InstitutionInviteButton institution={institution} />}
      <InstitutionBackgroundImageUploadDialog
        institution={institution}
        open={backgroundUploadOpen}
        onClose={() => setBackgroundUploadOpen(false)}
        onChange={(imgURL) => {
          setInstitution({
            ...institution,
            backgroundPictureURL: imgURL,
          });
        }}
      />
    </>
  );
}

export default InstitutionPage;
