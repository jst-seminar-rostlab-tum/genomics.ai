import React from 'react';
import Card from '@mui/material/Card';
import styles from './listCard.module.css';

function InstitutionCard({
  imageURL, title, description, nextToTitle, trailing,
}) {
  return (
    <Card sx={{ background: 'var(--card-background)' }}>
      <div style={{
        display: 'flex',
        flexDirection: 'row',
        alignItems: 'center',
        padding: '8px 16px',
      }}
      >
        <div className={styles.start}>
          {imageURL && (
          <img className={styles.cardImage} src={imageURL} alt="" />
          )}
          <div className={styles.text}>
            <div className={styles.titleRow}>
              <h3 className={styles.title}>{title}</h3>
              {nextToTitle || (<></>)}
            </div>
            <div className={styles.description}>
              {description && (
                <p>{description}</p>
              )}
            </div>
          </div>
        </div>
        {trailing && (
        <div className={styles.trailing}>
            {trailing}
        </div>
        )}
      </div>
    </Card>
  );
}

export default InstitutionCard;
